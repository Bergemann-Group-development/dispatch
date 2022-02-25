!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> A patch extension with pebbles.  There are so many patch properties that the
!> pebbles need, that it is much simpler to use an extension of an mhd_t, rather
!> than to pass (or/and recompute) those properties.  The key flag is do_pebbles,
!> which modifies slightly the variable allocations in mhd_mod.f90, so that
!> gas pressure and gas velocities are permanently allocated, and are available
!> to the pebble procedures.  To simplify neighbor access, the patch also has
!> an anonymous pointer to its primary list node, which has the nbor head of the
!> neighbor list.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pebble_patch_mod
  USE io_mod
  USE omp_lock_mod
  USE trace_mod
  USE bits_mod
  USE omp_mod
  USE omp_timer_mod
  USE link_mod
  USE kinds_mod
  USE task_mod
  USE patch_mod
  USE gpatch_mod
  USE mhd_mod
  !USE rmhd_mod
  USE planet_mod
  USE atmosphere_mod
  USE scaling_mod
  USE units_mod
  USE random_mod
  USE shared_mod
  USE pebble_mod
  USE pebbles_mod
  USE pebble_trace_mod
  USE mesh_mod
  USE solver_mod
  USE index_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! A patch with pebbles
  !-----------------------------------------------------------------------------
  type, public, extends(solver_t):: pebble_patch_t
    class(pebbles_t), pointer:: pebbles=>null()
    class(pebbles_t), pointer:: not_done=>null()
    procedure(kdk_drag), pointer:: move=>null() ! particle mover
    real:: pebble_rate                          ! pebble arrival rate
    real:: vz_pebble                            ! vertical in-fall speed
    logical:: only_pebbles                      ! freeze the MHD conditions (for tests)
  contains
    procedure init                              ! initialize the position and size
    procedure volume_add                        ! add particles throughout
    procedure boundary_add                      ! add particles at boundaries
    procedure drag                              ! compute gas drag
    procedure kdk_drag                          ! kick-drift-kick with drag
    procedure scatter                           ! scatter to not_done
    procedure gather                            ! gather from not_done
    procedure update_pebbles                    ! advance until self%time
    procedure info                              ! print info to stdout
    procedure index_r4                          ! compute interpolation index
    procedure inside_r4                         ! is it inside?
    procedure interpolate_r4                    ! interpolation
    procedure output                            ! pebble output
    procedure input                             ! pebble input
    procedure outside_sph
  end type
  logical, save:: debug = .false.
  integer, save:: id_dbg = 0
  integer, save:: verbose = 0
  integer, save:: n_vol=0, n_bdry=0
  integer, save:: track_level = 99
  real, save:: xmin=1e10, xmax=-1e10, z_max=1e9
  logical, save:: detailed_timer=.false.
  logical, save:: trace_any = .false.
  logical, save:: force_restart = .false.
  logical, save:: do_split=.true.
CONTAINS

!===============================================================================
!> Create a number of particles proportional at the left hand side wall
!===============================================================================
SUBROUTINE init (self)
  class(pebble_patch_t):: self
  !.............................................................................
  character(len=8), save:: type='kdk_drag'
  real, save:: density_yr_cgs = 1e-12, solid_density_cgs=3., size(3), pebble_size_cm=1.
  real, save:: pebble_rate_cgs = 0.0    ! 1e2 = about one pebble per 10 second (3e6 per year)
  real, save:: vz_pebble_cgs = 3e6      ! about 30 km/s vertically
  real, save:: size_range = 3.0         ! random number range in power of 10
  real, save:: ppc = 0.10               ! probability per cell of making a pebble
  real, save:: omega_yr(3)=0.0          ! omega in units where 1 = 1/yr
  integer, save:: trace_every=0         ! trace pebble frequency
  integer, save:: it_max=1000000, trace_every_min
  logical, save:: do_pebbles=.true., only_pebbles=.false., do_trace_output
  integer:: iostat, n(3)
  logical:: ok = .false.
  logical, save:: first_time=.true.
  real, dimension(:,:,:,:), pointer:: U
  namelist /pebble_params/ do_pebbles, type, density_yr_cgs, pebble_size_cm, z_max, &
    solid_density_cgs, size, size_range, it_max, pebble_rate_cgs, vz_pebble_cgs, &
    id_dbg, only_pebbles, ppc, omega_yr, trace_every, detailed_timer, do_trace_output, &
    trace_any, force_restart, do_split, track_level, verbose
  !-----------------------------------------------------------------------------
  ! Make sure to initialize, but only once
  !-----------------------------------------------------------------------------
  call trace%begin ('pebble_patch_t%init')
  allocate (self%pebbles,self%not_done)
  call self%pebbles%init (self%id)
  call self%not_done%init (self%id)
  !-----------------------------------------------------------------------------
  ! Read parameters, but only from one thread
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    call self%random%init
    self%random%seed = 11
    do_pebbles = .true.
    first_time = .false.
    do_trace_output = io%do_output
    size = self%size
    if (verbose>1) &
      print*, 'thread', omp%thread, ' reading namelist'
    rewind (io%input)
    read (io%input, pebble_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, pebble_params)
    end if
    write(io%namelist, pebble_params)
    if (do_pebbles) then
      n = self%n
      n(3) = min(n(3), ceiling(z_max*2./self%mesh(3)%d))
      trace_every_min = ceiling( &
        shared%n_patch_total*real(product(n))*ppc/real(pebble_trace%nbuffer))
      print *, 'trace_every_min =', trace_every_min, shared%n_patch_total
      if (trace_every < trace_every_min) then
        trace_every = trace_every_min
        print '("WARNING: trace_every has been increased to", i5)', trace_every
      end if
    end if
    call pebble_trace%init (trace_every)
  end if
  call io%lock%unset
  if (.not. do_pebbles) then
    call trace%end ()
    return
  end if
  self%omega = omega_yr*2.*cgs%pi*scaling%t/cgs%yr
  planet%omega = self%omega
  self%do_pebbles = do_pebbles
  self%only_pebbles = only_pebbles
  self%pebble_rate = pebble_rate_cgs*scaling%t
  self%vz_pebble = vz_pebble_cgs/scaling%u
  !-----------------------------------------------------------------------------
  ! Compute velocity and pressure
  !-----------------------------------------------------------------------------
  associate (mx=>self%mesh(1)%gn,my=>self%mesh(2)%gn,mz=>self%mesh(3)%gn,mt=>self%nt)
    self%density => self%mem(:,:,:,self%idx%d,:,1)
    allocate (self%pressure(mx,my,mz,mt))
    allocate (self%vgas(mx,my,mz,3,mt))
    U => self%vgas(:,:,:,:,self%it)
    call self%p2U (U, self%it)
    !print *,'pebble_patch_t%init: max(U) =',maxval(abs(U))
  end associate
  !-----------------------------------------------------------------------------
  ! Compute drag coefficient; the inverse stopping time for a 1 cm particle in
  ! a gas density of 1 g/cm**3, in the sub-sonic Epstein regime, specified in
  ! terms of the cgs gas density that gives one year stopping time.
  !-----------------------------------------------------------------------------
  self%pebbles%n_max = 3*ppc*product(self%n)
  self%pebbles%size_range = size_range
  self%pebbles%ppc = ppc
  self%pebbles%id_dbg = id_dbg
  self%pebbles%size = pebble_size_cm/scaling%l
  self%pebbles%it_max = it_max
  self%pebbles%density = solid_density_cgs/scaling%d
  self%pebbles%mass = self%pebbles%density*self%pebbles%size**3
  if (density_yr_cgs > 0.0) then
    self%pebbles%drag_coeff = scaling%d/(scaling%l*density_yr_cgs*cgs%yr/scaling%t)
  else
    self%pebbles%drag_coeff = 0.0
  end if
  pebble_trace%every = trace_every
  pebble_trace%do_trace_output = do_trace_output
  !-----------------------------------------------------------------------------
  ! Select mover
  !-----------------------------------------------------------------------------
  select case (type)
  case ('kdk_drag')
    self%move => kdk_drag
  end select
  !-----------------------------------------------------------------------------
  ! Try reading a restart snapshot and initialize new pebbles if that fails
  !-----------------------------------------------------------------------------
  if (.not.force_restart) call self%input (ok)
  if (ok) then
    ! ...
  else
    call self%boundary_add
    call self%volume_add
  end if
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Create a number of particles distributed randomly, with the velocity initially
!> equal to the gas velocity.  To be able to use single precision particle
!> coordinates, they are measured relative to the position of the planet.
!===============================================================================
SUBROUTINE volume_add (self, n)
  class(pebble_patch_t):: self
  integer, intent(in), optional:: n
  !.............................................................................
  class(pebbles_t), pointer:: pebbles
  class(pebble_t), pointer:: pebble
  integer:: ix, iy, iz, ippc, iter, warn, l(3), u(3), i, ncell, npart
  real:: w(3), sign, ppc, average_density, md(3)
  real(kind=8):: p(3)
  integer, save:: nprint=10
  class(mesh_t), pointer:: m
  !-----------------------------------------------------------------------------
  if (.not.self%do_pebbles .or. self%is_set(bits%not_leaf)) return
  do i=1,3
    md(i)=max((self%mesh(i)%lf+self%position(i))**2,(self%mesh(i)%uf+self%position(i))**2)
  end do
  !if (sum(md) >= planet%bondi_radius**2) then
    m => self%mesh(3)
    !---------------------------------------------------------------------------
    ! If the patch lies entirely above or below +- z_max, skip it
    !---------------------------------------------------------------------------
    if ((m%p-m%s*0.5) > z_max .or. (m%p+m%s*0.5) < -z_max) return
    call trace%begin ('pebbles_t%volume_add', 3)
    pebbles => self%pebbles
    if (verbose>1 .or. (verbose>0.and.nprint>0)) then
      nprint = nprint-1
      print *, ' ADDING RANDOM PEBBLES TO PATCH', self%id
      print 1, '     patch position:', self%position - planet%position
      print 1, '         patch size:', self%size
      print 1, '   pebble size (cm):', self%pebbles%size*scaling%l
      print 1, 'solid density (cgs):', self%pebbles%density*scaling%d
      print 1, '         mass (cgs):', self%pebbles%mass*scaling%m
      print 1, '   model size (cgs):', self%size*scaling%l
      print 1, ' gas velocity (cgs):', self%vgas(3,3,3,1:3,self%it)*scaling%u
      print 1, '      pebbles added:', self%pebble_rate*self%dtime, self%pebble_rate, self%dtime
      1 format(a,1p,3e11.2)
    end if
    l = self%mesh%li
    u = self%mesh%ui
    average_density = sum(self%density(l(1):u(1),l(2):u(2),l(3):u(3),self%it))/product(u-l+1)
    ncell = 0
    npart = 0
    do iz=self%mesh(3)%li,self%mesh(3)%ui
      m => self%mesh(3)
      p(3) = m%p + m%r(iz) - planet%position(3)
      !-------------------------------------------------------------------------
      ! If the cell lies entirely below or above +- z_max, skip it
      !-------------------------------------------------------------------------
      if (p(3)-0.5*m%d > +z_max) cycle
      if (p(3)+0.5*m%d < -z_max) cycle
      do iy=self%mesh(2)%li,self%mesh(2)%ui
        m => self%mesh(2)
        p(2) = m%p + m%r(iy) - planet%position(1)
        do ix=self%mesh(1)%li,self%mesh(1)%ui
          !---------------------------------------------------------------------
          ! Make sure the location is outside the planet radius
          !---------------------------------------------------------------------
          !if (self%outside_sph([ix,iy,iz])) then
          m => self%mesh(1)
          p(1) = m%p + m%r(ix) - planet%position(1)
          if (sum(p**2) < planet%radius**2) cycle
          ncell = ncell+1
          !---------------------------------------------------------------------
          ! Probability ia uniform, but weight is proportional to local density
          !---------------------------------------------------------------------
          ppc = self%pebbles%ppc
          ippc = ppc
          ppc = ppc-ippc
          do iter=1,ippc+1
            if (iter<=ippc .or. (ppc>0.0 .and. self%random%ran1()<ppc)) then
              allocate (pebble)
              npart = npart+1
              !$omp atomic
              n_vol = n_vol+1
              call pebble%init
              if (pebble_trace%every > 0) then
              if (mod(pebble%id,pebble_trace%every)==0) &
                pebble%status = ior(pebble%status,bits%trace+bits%has_BC)
              end if
              pebble%t = self%time
              pebble%s = self%pebbles%size*10.0**(-self%pebbles%size_range*self%random%ran1())
              pebble%m = pebble%s**3*self%pebbles%density
              !-----------------------------------------------------------------
              ! Initially place the particles uniformly in the cell
              !-----------------------------------------------------------------
              pebble%p(1) = p(1) + m%d*(0.5-self%random%ran1())
              pebble%p(2) = p(2) + m%d*(0.5-self%random%ran1())
              pebble%p(3) = p(3) + m%d*(0.5-self%random%ran1())
              !-----------------------------------------------------------------
              ! If the cell is higher than 2*z_max, scale down the z-range.
              ! FIXME: This may give slightly wrong results unless the affected
              ! cell is centered around z=0.0
              !-----------------------------------------------------------------
              if (self%mesh(3)%d > z_max*2.0) then
                pebble%p(3) = pebble%p(3)*z_max*2.0/self%mesh(3)%d
              end if
              if (sqrt(sum(pebble%p**2)) < 1e-10) print *, pebble%id, pebble%p
              !-----------------------------------------------------------------------
              ! Avoid placement outside (could happen for no-no-mans-land mesh)
              !-----------------------------------------------------------------------
              if (any(pebble%p < self%mesh%p+self%mesh%lf .or. &
                      pebble%p > self%mesh%p+self%mesh%uf)) then
                !print *, pebble%p, real(self%mesh%p+self%mesh%lf), real(self%mesh%p+self%mesh%uf)
                deallocate(pebble)
                cycle
              end if
              !-----------------------------------------------------------------------
              ! Single precision position measured relative to planet position, with
              ! initial speed = gas speed, and weight proportional to local density
              !-----------------------------------------------------------------------
              call interpolate_vg (self, pebble, warn)
              if (self%pebbles%drag_coeff > 0.0) then
                pebble%v = pebble%vg
              else
                pebble%v = 0.0
              end if
              pebble%d = self%mem(ix,iy,iz,self%idx%d,self%it,1)
              pebble%w = pebble%d*product(self%ds)/self%pebbles%ppc
              if (pebble_trace%every > 0) then
                if (trace_any.and.mod(pebble%id,pebble_trace%every)==0) then
                  pebble%status = ior(pebble%status,bits%trace)               ! trace some
                  !$omp atomic
                  shared%n_initial = shared%n_initial+1
                else
                  if (mod(pebble%id,pebble_trace%every)==0 .and. &
                abs(pebble%p(1)) < 2.0*planet%hill_radius .and. &
                abs(pebble%p(2)) < 2.0*planet%hill_radius) then
                    pebble%status = ior(pebble%status,bits%trace)               ! trace some
                    !$omp atomic
                    shared%n_initial = shared%n_initial+1
                  end if
                end if
              end if
              pebble%t = self%time                                            ! time
              call self%pebbles%append (pebble)
            end if
          end do
        end do
      end do
    end do
    if (verbose>0 .and. ncell>0) &
      print *,self%id, 'ADD: ncell, npart =', ncell, npart, self%pebbles%n, &
              real(self%size(3)/2.), real(self%position(3))
    if (verbose>2) &
      print '(a,i7,i3,2i6,1p,2e14.6)', 'pebble_patch%volume_add: id, level, n_initial, n_add =', &
        self%id, self%level, shared%n_initial, shared%n_add, xmin, xmax
    if (verbose>1) &
      print *, 'volume_add:', self%id, self%position(3), self%pebbles%n
    call self%pebbles%check_list ('after volume_add')
    call trace%end
  !end if
END SUBROUTINE volume_add

LOGICAL FUNCTION outside_sph(self, cl)
  class(pebble_patch_t) :: self
  integer :: cl(3)
  ! ....................................................................
  class(mesh_t), pointer:: m
  real(8) :: pos(3), dst
  integer :: i
  ! --------------------------------------------------------------------
  outside_sph = .true.
  dst = 0d0
  do i=1,3
    m => self%mesh(i)
    pos(i) = m%r(cl(i))+self%position(i)
    dst = dst + (pos(i) - planet%position(i))**2
  end do
  if (dst < planet%bondi_radius**2) outside_sph = .false.
END FUNCTION outside_sph

!===============================================================================
!> Create a number of particles proportional at the top and bottom walls,
!> with initially a velocity interpolated from the argument velocity field.
!> Depending on the call, either add a fixedd number, or a number proportional
!> to the time step (and hence very much larger in the large, cold, outer
!> patches.  Building the rate into the weight factor may be a better strategy.
!===============================================================================
SUBROUTINE boundary_add (self, n)
  class(pebble_patch_t):: self
  integer, intent(in), optional:: n
  !.............................................................................
  class(pebbles_t), pointer:: pebbles
  class(pebble_t), pointer:: pebble
  integer:: i, j, k, ii(3), n_local
  real:: w(3), sign, x0
  integer, save:: nprint=10
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.self%do_pebbles .or. self%pebble_rate==0.0) return
  if (self%level>1 .or. self%is_clear(bits%top_bottom)) return
  call trace%begin ('pebbles_t%boundary_add', 3, itimer=itimer)
  pebbles => self%pebbles
  if (present(n)) then
    n_local = n
  else
    n_local = self%pebble_rate*self%dtime
  end if
  if (verbose>1 .or. (verbose>0.and.nprint>0)) then
    nprint = nprint-1
    print *, ' ADDING TOP/BOTTOM PEBBLES TO PATCH', self%id
    print 1, '     patch position:', self%position - planet%position
    print 1, '         patch size:', self%size
    print 1, '   pebble size (cm):', self%pebbles%size*scaling%l
    print 1, 'solid density (cgs):', self%pebbles%density*scaling%d
    print 1, '         mass (cgs):', self%pebbles%mass*scaling%m
    print 1, '   model size (cgs):', self%size*scaling%l
    print 1, ' gas velocity (cgs):', self%vgas(3,3,3,1:3,self%it)*scaling%u
    print 1, '      pebbles added:', self%pebble_rate*self%dtime, self%pebble_rate, self%dtime
    1 format(a,1p,3e11.2)
  end if
  self%pebbles%z_min = 1e30
  do i=1,n_local
    allocate (pebble)
    call pebble%init
    pebble%t = self%time
    !$omp atomic
    n_vol = n_vol+1
    pebble%s = self%pebbles%size*10.**(-self%pebbles%size_range*self%random%ran1())
    pebble%m = pebble%s**3*pebble%d
    sign = merge(1.0,-1.0,self%position(3)>0.0)
    pebble%p(1) =   (self%position(1) - planet%position(1)) + self%size(1)*(0.5-self%random%ran1())
    pebble%p(2) =   (self%position(2) - planet%position(2)) + self%size(2)*(0.5-self%random%ran1())
    pebble%p(3) =   (self%position(3) - planet%position(3)) + self%size(3)*0.5*sign
    pebble%v(1) = - self%omega(3)*planet%orbit_radius
    pebble%v(2) =   0.0
    pebble%v(3) = - self%vz_pebble*sign
    ii = self%index_r4 (pebble%p, w, self%mesh%gn)
    pebble%d = self%interpolate_r4 (pebble%t, self%density , ii, w)
    pebble%w = product(self%size)/n_local                               ! probability density
    pebble%w = pebble%w*pebble%d                                        ! normalize by gas density
    self%pebbles%z_min = min(self%pebbles%z_min,abs(pebble%p(3)))
    !---------------------------------------------------------------------------
    ! If the particle has a chance to reach the planet it is a candidate for
    ! tracing -- select every pebble_trace%every particle for tracing
    !---------------------------------------------------------------------------
    x0 = pebble%p(1) !- pebble%v(1)*pebble%p(2)/pebble%v(2)
    if (trace_any.and.mod(pebble%id,pebble_trace%every)==0) then
      pebble%status = ior(pebble%status,bits%trace)
      !$omp atomic
      shared%n_add = shared%n_add+1
      xmin = min(xmin,x0)
      xmax = max(xmax,x0)
    else if (x0 > -3.0*planet%hill_radius .and. &
             x0 <  3.0*planet%hill_radius .and. &
             abs(pebble%p(2)) < 5.0*planet%hill_radius .and. &
             pebble_trace%every > 0) then
      if (mod(pebble%id,pebble_trace%every)==0) then
        pebble%status = ior(pebble%status,bits%trace)
        !$omp atomic
        shared%n_add = shared%n_add+1
        xmin = min(xmin,x0)
        xmax = max(xmax,x0)
      end if
    end if
    if (debug.or.pebble%id==id_dbg) print*,'dbg: add', i, n_local, pebble%id
    call self%pebbles%append (pebble)
  end do
  call self%pebbles%check_list ('after add')
  call trace%end (itimer)
END SUBROUTINE boundary_add

!===============================================================================
!> Look up gas velocity components, given positions, and compute the drag rate
!> (inverse stopping time)
!===============================================================================
SUBROUTINE drag (self, pebble, warning)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: pebble
  !.............................................................................
  real:: w(3), pos(3), d, p, temp, gamma, v(3), speed2, sound2
  integer:: ii(3), j, k
  integer, save:: nprint=10
  integer:: w1, w2, warning
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ii  = self%index_r4 (pebble%p, w, self%mesh%gn)
  d   = self%interpolate_r4 (pebble%t, self%density , ii, w)
  p   = self%interpolate_r4 (pebble%t, self%pressure, ii, w)
  call interpolate_vg (self, pebble, w2)
  if (detailed_timer) call trace%begin ('pebble_patch_t%drag', itimer=itimer)
  w1  = count(w<0.0 .or. w>1.0)
  warning = w1+w2
  if (pebble%id==id_dbg) print '(i8,2x,a,i6,1p,3(3x,3e10.2))', &
    id_dbg, 'interpolate_r4: id, p, vg, d', self%id, pebble%p, pebble%vg, d
  !
  speed2    = sum((pebble%v-pebble%vg)**2)
  temp      = p/d
  gamma     = merge(1.6666667,1.4,temp>2000.)
  sound2    = gamma*p/d
  pebble%d  = d
  pebble%c  = self%pebbles%drag_coeff*d/pebble%s*sqrt((speed2+sound2)/sound2)
  pebble%ad = (pebble%vg-pebble%v)*pebble%c
  if (warning>0 .and. verbose>2) then
    pos = pebble%p + real(self%position-planet%position,kind=4)
    print '("ww",2i1,1p,9(3e10.2,1x))', &
      w1,w2, &
      pos/self%ds, &
      (pos+2.*pebble%dt*pebble%v)/self%ds, &
      pebble%v*pebble%dt/self%ds, &
      pebble%ag, pebble%dt, pebble%c*pebble%dt
  end if
  if (verbose>0 .and. nprint>0) then
    nprint = nprint-1
    print *,'                  id:', pebble%id
    print *,'    pebble size (cm):', pebble%s*scaling%l
    print *,'  stopping time (yr):', 1./pebble%c*scaling%t/cgs%yr
    print *,'   gas density (cgs):', d*scaling%d
    print *,'  gas velocity (cgs):', pebble%vg*scaling%u
    print *,'      velocity (cgs):', pebble%v*scaling%u
    print *,'     temperature (K):', temp*scaling%temp
    print *,'   sound speed (cgs):', sqrt(sound2)*scaling%u
    print *,'relative speed (cgs):', sqrt(speed2)*scaling%u
    print *,'   drag acceleration:', pebble%ad*scaling%l/scaling%t**2
    print *,'gravity acceleration:', pebble%ag*scaling%l/scaling%t**2
    print *,'       time step (s):', pebble%dt*scaling%t,pebble%dt
    print *,'           stiffness:', pebble%dt*pebble%c
    print *,'     cells traversed:', pebble%dt*pebble%v/self%ds
  end if
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE drag

!===============================================================================
!> Move a pebbles forward one step, with drag
!===============================================================================
SUBROUTINE kdk_drag (self, p)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: p
  !.............................................................................
  real:: rp(3), ar, ds
  real(8) :: rad, sp, p1(3), ag(3)
  integer, save:: nprint=20
  integer:: warn
  integer, save:: itimer=0
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Get acceleration from planet and Coriolis force from star, ag(t)
  !---------------------------------------------------------------------------
  if (iand(p%status,bits%frozen) /= 0) return
  if (detailed_timer) call trace%begin ('pebble_patch_t%kdk_drag', itimer=itimer)
  if (planet%r_smooth>0) then
    rad = sqrt(sum(p%p**2))
    sp = p%p(1); p%ag(1) = planet%newton(rad,sp)
    sp = p%p(2); p%ag(2) = planet%newton(rad,sp)
    sp = p%p(3); p%ag(3) = planet%newton(rad,sp)
  else
    rp = - p%p
    ar = 1./max(sqrt(sum(rp**2)),planet%radius*0.5)
    p%ag = scaling%grav*planet%mass*rp*ar**3
  end if
  p%ag(1) = p%ag(1) + 2.*self%omega(3)*p%v(2) + self%omega(3)**2*p%p(1)*3.0
  p%ag(2) = p%ag(2) - 2.*self%omega(3)*p%v(1)
  p%ag(3) = p%ag(3)                           - self%omega(3)**2*p%p(3)
  !rp = self%pebbles%position - p%p
  !ar = 1./sqrt(sum(rp**2))
  !p%ag = p%ag + scaling%grav*star%mass*rp*ar**3
  ag = p%ag
  p1 = p%p
  p%dt = sqrt(sum(p1**2)/sqrt(sum(ag**2)))*self%pebbles%courant
  p%dt = min(p%dt,0.1*io%out_time,minval(self%ds/max(1e-20,abs(p%v))))
  p%dt = min(p%dt,0.51_8*(self%time-p%t))
  !
  if (p%id==self%pebbles%id_dbg) then
    ds = sqrt(sum((p%dt*p%v/self%ds)**2))
    print'("kdkd0: it, max(vg), p%t, t(iit) =",i6,1p,e14.5,2x,7e14.5)', &
      self%id,maxval(abs(self%vgas(:,:,:,:,self%it))),p%t,self%t(self%iit)
    print *,'kdkd: p,v,vg,ad,cdt,dt,t,ds '
    print'("kdkd1",1p,7(3e10.2,1x))', &
      p%p, p%v, p%vg, p%ad, p%c*p%dt, p%dt, p%t, ds
  end if
  !
  if (verbose>0 .and. p%id==1 .and. nprint>0) then
    nprint = nprint-1
    print'("kdkd:p,v,ag,ds,t,dt",1p,6(3e10.2,2x))', &
      p%p, p%v, p%ag, self%ds, p%t, p%dt
  end if
  !---------------------------------------------------------------------------
  ! Get the drag acceleration, at t+dt, then back up position. This requires
  ! access to patch data at time=p%t+p%dt, but not access to data at time=p%t.
  ! Hence we can wait with moving particle until (at least) the time when
  ! p%t < self%t(self%iit(1))
  !---------------------------------------------------------------------------
  p%p = p%p + p%dt*p%v
  p%t = p%t + p%dt
  call self%drag (p,warn)
  if (p%id==self%pebbles%id_dbg) then
    ds = sqrt(sum((p%dt*p%v/self%ds)**2))
    print'("kdkd2",1p,7(3e10.2,1x))', p%p, p%v, p%vg, p%ad, p%c*p%dt, p%dt, p%t, ds
  end if
  if (warn>0 .and. verbose>3) then
    print *,'  dt',p%dt,1,verbose
    print *,' dt1',self%pebbles%courant*sqrt(sum(p%v**2)/sum(p%ag**2))
    print *,' dt2',self%pebbles%courant*sqrt(sum(p%p**2)/sqrt(sum(p%ag**2)))
    print *,' dt3',minval(self%ds)/max(1e-20,abs(p%v)),io%out_time
    print'(1x,a,1p,9(3e11.2,2x))', 'dt3', self%ds, abs(p%v), self%ds/max(1e-20,abs(p%v))
    print *,'  ds',self%ds
    print *,'v*dt',p%v*p%dt
  end if
  p%p = p%p - p%dt*p%v
  !---------------------------------------------------------------------------
  ! First kick, with drag, v(t) -> v(t+dt)
  !---------------------------------------------------------------------------
  p%v = p%v + p%dt*(p%ag+p%ad)/(1.+p%c*p%dt)
  !---------------------------------------------------------------------------
  ! Drift, r -> r+2*dt*v(t+dt*a), with v -> vg(t+dt) if drag dominates
  !---------------------------------------------------------------------------
  p%p = p%p + 2.*p%dt*p%v
  p%t = p%t + p%dt
  !---------------------------------------------------------------------------
  ! Get acceleration from planet and Coriolis force from star, at p(t+2*dt)
  !---------------------------------------------------------------------------
  rp = - p%p
  ar = 1./max(sqrt(sum(rp**2)),planet%radius*0.5)
  p%ag = scaling%grav*planet%mass*rp*ar**3
  !rp = self%pebbles%position - p%p
  !ar = 1./sqrt(sum(rp**2))
  !p%ag = p%ag + scaling%grav*star%mass*rp*ar**3
  if (first_time) print *, 'kdk_drag: omega =', self%omega
  p%ag(1) = p%ag(1) + 2.*self%omega(3)*p%v(2) + self%omega(3)**2*p%p(1)*3.0
  p%ag(2) = p%ag(2) - 2.*self%omega(3)*p%v(1)
  p%ag(3) = p%ag(3)                           - self%omega(3)**2*p%p(3)
  !---------------------------------------------------------------------------
  ! Get the drag acceleration, at t+2*dt
  !---------------------------------------------------------------------------
  call self%drag (p,warn)
  if (warn>0 .and. verbose>3) then
    print *,'  dt',p%dt,2
    print *,' dt1',self%pebbles%courant*sqrt(sum(p%v**2)/sum(p%ag**2))
    print *,' dt2',minval(self%ds/max(1e-20,abs(p%v)))
    print *,' dt3',min(p%dt,minval(self%ds/max(1e-20,abs(p%v))))
    print'(1x,a,1p,9(3e11.2,2x))', 'dt3', self%ds, abs(p%v), self%ds/max(1e-20,abs(p%v))
    print *,'  ds',self%ds
    print *,'v*dt',p%v*p%dt
  end if
  !---------------------------------------------------------------------------
  ! Second kick, with drag, moving v to v(t+2*dt)
  !---------------------------------------------------------------------------
  if (p%id==self%pebbles%id_dbg) then
    ds = sqrt(sum((p%dt*p%v/self%ds)**2))
    print'("kdkd3",1p,7(3e10.2,1x))', p%p, p%v, p%vg, p%ad, p%c*p%dt, p%dt, p%t, ds
  end if
  p%v = p%v + p%dt*(p%ag+p%ad)/(1.+p%c*p%dt)
  if (p%id==self%pebbles%id_dbg) then
    nprint = nprint-1
    ds = sqrt(sum((p%dt*p%v/self%ds)**2))
    print'("kdkd4",1p,7(3e10.2,1x))', p%p, p%v, p%vg, p%ad, p%c*p%dt, p%dt, p%t, ds
  end if
  if (p%id == id_dbg) &
    print*,'dbg: kdk_drag:', p%p
  if (detailed_timer) call trace%end (itimer)
  if (iand(p%status,bits%trace)/=0) then
    call pebble_trace%trace_pebble (p, self%id, self%level)
  end if
  !$omp atomic write
  first_time = .false.
END SUBROUTINE kdk_drag

!===============================================================================
!> Starting from a collection of particles belonging to this patch, kept in the
!> self%pebbles object, process them one by one.  Those that are not yet behind
!> the patch time slices are left alone, while other ones are moved one step,
!> after which it is checked if they should be exported to another patch.
!===============================================================================
SUBROUTINE scatter (self, nbors, not_done)
  class(pebble_patch_t):: self
  class(nbor_t), pointer:: nbors
  !.............................................................................
  class(pebbles_t), pointer:: not_done
  class(pebble_t), pointer:: pebble, next
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task
  integer:: n, m
  real(8):: pos(3)
  real:: w_norm
  real, pointer:: d(:,:,:)
  type(index_t):: idx
  integer, save:: itimer=0
  real, dimension(3):: p_save, v_save
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('pebble_patch_t%scatter', itimer=itimer)
  self%pebbles%dt_min = 1e30
  self%pebbles%dt_max = 0.0
  self%pebbles%z_min = 1e30
  n = 0
  nullify(not_done%head)
  nullify(not_done%tail)
  pebble => self%pebbles%head
    if (debug.or.pebble%id==id_dbg) print *,'scatter: patch pos, size:', pebble%id, self%position, self%size
  do while (associated(pebble) .and. self%pebbles%n>0)
    !---------------------------------------------------------------------------
    ! Lock access to the pebble list while getting the next pointer, since some
    ! other patch may be trying exporting to us at the same time.
    !---------------------------------------------------------------------------
    call self%pebbles%lock%set ('scatter next')
    next => pebble%next
    !$omp flush
    call self%pebbles%lock%unset ('scatter next')
    !---------------------------------------------------------------------------
    ! If the pebble is already ahead in time, or frozen, just skip it.  If we
    ! use >= here, particles will not move until all time slots have advanced.
    !---------------------------------------------------------------------------
    if (pebble%t>=self%t(self%iit(1)) .or. iand(pebble%status,bits%frozen) /= 0) goto 9
    n = n+1
    if (verbose>0) then
      p_save = pebble%p
      v_save = pebble%v
      call self%move (pebble)
      call check_nan
    else
      call self%move (pebble)
    end if
    self%pebbles%z_min = min(self%pebbles%z_min,abs(pebble%p(3)))
    self%pebbles%dt_min = min(pebble%dt,self%pebbles%dt_min)
    self%pebbles%dt_max = max(pebble%dt,self%pebbles%dt_max)
    if (debug.or.pebble%id==id_dbg) print *, pebble%id, 'scatter: pebble id:', self%is_set(bits%not_leaf), pebble%p
    !---------------------------------------------------------------------------
    ! Check if it should be exported to a neighbor, either because it is
    ! outside, or because there is a higher resolution child patch.
    !---------------------------------------------------------------------------
    if (.not. self%inside_r4 (pebble%p) .or. self%is_set(bits%not_leaf)) then
      if (debug.or.pebble%id==id_dbg) print *, pebble%id, 'scatter: not inside its patch', self%id
      if (associated(nbors)) then
        m = 0
        nbor => nbors
        do while (associated(nbor))
          task => nbor%task
          select type (task)
          class is (pebble_patch_t)
            if (task%inside_r4(pebble%p)) then
              if (.not.(self%is_set(bits%not_leaf)) .or. task%level>self%level) then
                if (debug.or.pebble%id==id_dbg) &
                  print*, pebble%id, 'dbg: exporting to other task', task%id, task%pebbles%id, task%pebbles%n
                !---------------------------------------------------------------
                ! Periodic wrapping, in case the new patch is a periodic nbor
                !---------------------------------------------------------------
                pos = pebble%p + planet%position
                where (self%periodic)
                  pos = modulo(pos+self%box*0.5_8,self%box) - self%box*0.5_8
                end where
                pebble%p = pos - planet%position
                if ((verbose > 0) .and. any(abs(self%position-task%position) > self%box)) then
                  print *, 'INFO: exporting across periodic boundary', self%id, task%id, pebble%id
                end if
                if (do_split) then
                  d => task%mem(:,:,:,idx%d,task%it,1)
                  w_norm = task%faver(d)*product(task%ds)/self%pebbles%ppc
                  call self%pebbles%export (pebble, task%pebbles, w_norm, task%ds)
                else
                  call self%pebbles%export (pebble, task%pebbles)
                end if
                goto 9
              end if
            end if
          end select
          nbor => nbor%next
          m = m+1
          if (m > 200) then
            print*, self%id, 'ERROR: stuck in nbor loop'
            exit
          end if
        end do
        !-----------------------------------------------------------------------
        ! If not inside, and not inside a neighbor, remove it
        !-----------------------------------------------------------------------
        if (debug.or.pebble%id==id_dbg) &
          print *,pebble%id,'scatter: not inside nbor: remove'
        call pebble_trace%check (pebble)
        call self%pebbles%remove (pebble)
        self%pebbles%n_remove = self%pebbles%n_remove + 1
        if (verbose>0) &
          call verbose_remove
        goto 9
      else
        !-----------------------------------------------------------------------
        ! No neighbors???
        !-----------------------------------------------------------------------
        if (debug.or.pebble%id==id_dbg) &
          print *,'WARNING: pebble_patch without neighbors', self%id
        if (abs(pebble%p(1)) > self%size(1)*0.5) &
          self%pebbles%n_outside(1) = self%pebbles%n_outside(1) + 1
        if (abs(pebble%p(2)) > self%size(2)*0.5) &
          self%pebbles%n_outside(2) = self%pebbles%n_outside(2) + 1
        if (abs(pebble%p(3)) > self%size(3)*0.5) &
          self%pebbles%n_outside(3) = self%pebbles%n_outside(3) + 1
        if (debug.or.pebble%id==id_dbg) print *, pebble%id, 'scatter: no nbor: remove'
        call pebble_trace%check (pebble)
        call self%pebbles%remove (pebble)
        self%pebbles%n_remove = self%pebbles%n_remove + 1
        if (verbose>0) &
          call verbose_remove
        goto 9
      end if
    else if (planet%inside(pebble%p)) then
      !-------------------------------------------------------------------------
      ! The pebble hit the planet; accreted
      !-------------------------------------------------------------------------
      if (debug.or.pebble%id==id_dbg) &
        print *, self%id, 'scatter: inside planet: accrete', pebble%id
      pebble%status = ior(pebble%status,bits%frozen)
      call pebble_trace%check (pebble)
      call self%pebbles%remove (pebble)
      self%pebbles%n_accrete = self%pebbles%n_accrete + 1
      goto 9
    else if (pebble%t < self%time) then
      !-------------------------------------------------------------------------
      ! Finally, check if the pebble is inside a neighbor with higher level
      !-------------------------------------------------------------------------
      if (associated(nbors)) then
        m = 0
        nbor => nbors
        do while (associated(nbor))
          task => nbor%task
          if (task%level==self%level+1) then
            select type (task)
            class is (pebble_patch_t)
              if (task%inside_r4(pebble%p)) then
                if (debug.or.pebble%id==id_dbg) &
                  print*, pebble%id, 'dbg: exporting to other task', task%id, task%pebbles%id, task%pebbles%n
                if (do_split) then
                  d => task%mem(:,:,:,idx%d,task%it,1)
                  w_norm = task%faver(d)*product(task%ds)/self%pebbles%ppc
                  call self%pebbles%export (pebble, task%pebbles, w_norm, task%ds)
                else
                  call self%pebbles%export (pebble, task%pebbles)
                end if
                goto 9
              end if
            end select
          end if
          nbor => nbor%next
          m = m+1
          if (m > 200) then
            print*, self%id, 'ERROR: stuck in nbor loop'
            exit
          end if
        end do
      end if
      !-------------------------------------------------------------------------
      ! If it wasn't, export to the 'not_done' list
      !-------------------------------------------------------------------------
      if (debug.or.pebble%id==id_dbg) &
        print*, self%id, 'dbg: exporting to not_done', pebble%id, not_done%id, not_done%n
        call self%pebbles%export (pebble, not_done)
    else
      if (debug.or.pebble%id==id_dbg) print*, self%id, &
       'dbg: done, time, self%time =', pebble%id, pebble%t, self%time
    end if
  9 continue
    pebble => next
!    if (n > 50*self%pebbles%n) then
!      print*,'ERROR: stuck in scatter loop, moves, pebbles%n =', n, self%pebbles%n
!      exit
!    end if
  end do
  !if (debug) &
  if ((self%track .or. self%level>=track_level) .and. n>0) then
    print '(i5,1x,a,7i7,1p,2e12.3)', self%id, &
      'scatter: pebbles, accrete, remove, split, import, export, n', &
      self%pebbles%n, self%pebbles%n_accrete, self%pebbles%n_remove, self%pebbles%n_split, &
      self%pebbles%n_import, self%pebbles%n_export, n, self%pebbles%z_min
      self%pebbles%n_accrete = 0
      self%pebbles%n_import = 0
      self%pebbles%n_export = 0
      self%pebbles%n_remove = 0
      self%pebbles%n_frozen = 0
      self%pebbles%n_split  = 0
  end if
  call self%pebbles%check_list ('after scatter')
  !if (n > self%pebbles%n + &
  !        self%not_done%n + &
  !        self%pebbles%n_export + &
  !        self%pebbles%n_accrete + &
  !        self%pebbles%n_remove) then
  !  print*, self%id, ' ERROR: wrong pebble count, n, pebbles%n not_done%n =', &
  !   n, self%pebbles%n, self%not_done%n, self%pebbles%n_export, self%pebbles%n_accrete, self%pebbles%n_remove
  !end if
  if (verbose>1 .and. self%pebbles%n>1) then
    print*,'scatter: n, self%pebbles%n, not_done=', self%pebbles%n, n, self%not_done%n
    self%pebbles%n_update = self%pebbles%n_update + self%pebbles%n
  end if
  if (self%pebbles%n==0) then
    nullify(self%pebbles%head, self%pebbles%tail)
  end if
  if (detailed_timer) call trace%end (itimer)
contains
  subroutine verbose_remove
    logical:: inside
    nbor => nbors
    print *, 'verbose_remove', pebble%p, isnan(pebble%p)
    inside = self%inside_r4 (pebble%p, debug=.true.)
    do while (associated(nbor))
      task => nbor%task
      select type (task)
      class is (pebble_patch_t)
        inside = task%inside_r4 (pebble%p, debug=.true.)
      end select
      nbor => nbor%next
    end do
  end subroutine verbose_remove
  subroutine check_nan
    if (any(isnan(pebble%p))) then
      print *, 'NAN:', p_save, v_save
    end if
  end subroutine check_nan
END SUBROUTINE scatter

!===============================================================================
!===============================================================================
SUBROUTINE gather (self, not_done, nbors, verbose)
  class(pebble_patch_t):: self
  class(pebbles_t), pointer:: not_done
  class(nbor_t), pointer:: nbors
  !.............................................................................
  class(nbor_t), pointer:: nbor
  class(pebble_t), pointer:: pebble, next
  class(task_t), pointer:: task
  integer:: n, verbose
  real:: w_norm
  real, pointer:: d(:,:,:)
  type(index_t):: idx
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('pebble_patch_t%gather', itimer=itimer)
  not_done%dt_min = 1e30
  not_done%dt_max = 0.0
  not_done%t_min = 1e30
  n = 0
  pebble => not_done%head
  if (verbose>0 .or. debug.or.pebble%id==id_dbg) print*, self%id, &
    'gather: self%n, not_done%n =', self%pebbles%n, not_done%n
  !-----------------------------------------------------------------------------
  ! This loop moves each particle one step, and checks if has gone out of bounds
  !-----------------------------------------------------------------------------
  do while (associated(pebble) .and. not_done%n>0)
    next => pebble%next
    n = n+1
    call self%move (pebble)
    !self%pebbles%t_min  = min(pebble%t ,self%pebbles%t_min)
    !self%pebbles%dt_min = min(pebble%dt,self%pebbles%dt_min)
    !self%pebbles%dt_max = max(pebble%dt,self%pebbles%dt_max)
    if (.not. self%inside_r4 (pebble%p) .or. self%is_set(bits%not_leaf)) then
      !-------------------------------------------------------------------------
      ! Check if it should be exported to a neighbor, either because it is
      ! outside, or because there is a higher resolution child patch.
      !-------------------------------------------------------------------------
      nbor => nbors
      do while (associated(nbor))
        task => nbor%task
        select type (task)
        class is (pebble_patch_t)
          if (task%inside_r4(pebble%p)) then
            if (.not.(self%is_set(bits%not_leaf)) .or. task%level>self%level) then
              if (verbose>0 .or. debug.or.pebble%id==id_dbg) print*, self%id, &
                'gather: exporting to other task', pebble%id, task%id, task%pebbles%id, task%pebbles%n
              if (do_split) then
                d => task%mem(:,:,:,idx%d,task%it,1)
                w_norm = task%faver(d)*product(task%ds)/self%pebbles%ppc
                call not_done%export (pebble, task%pebbles, w_norm, task%ds)
              else
                call not_done%export (pebble, task%pebbles)
              end if
              not_done%n_export = not_done%n_export+1
              goto 9
            end if
          end if
        end select
        nbor => nbor%next
      end do
      call not_done%remove (pebble)
      not_done%n_remove = not_done%n_remove + 1
      if (verbose>0 .or. debug.or.pebble%id==id_dbg) print*, self%id, 'gather: no nbor: removed'
    else if (planet%inside(pebble%p)) then
      pebble%status = ior(pebble%status,bits%frozen)
      call pebble_trace%check (pebble)
      call not_done%remove (pebble)
      not_done%n_remove = not_done%n_remove + 1
      if (verbose>0 .or. debug.or.pebble%id==id_dbg) print *, self%id, 'gather: accreted'
    else if (pebble%t >= self%time) then
      call not_done%export (pebble, self%pebbles)
      not_done%n_export = not_done%n_export+1
      if (verbose>0 .or. debug.or.pebble%id==id_dbg) print*, self%id,&
       'gather: exporting to self', self%pebbles%id, self%pebbles%n, self%not_done%n
    end if
    if (verbose>0) print*, self%id, 'pebble%next', pebble%id, associated(pebble%next), pebble%t, pebble%dt, self%time, iand(pebble%status,bits%frozen) /= 0
  9 pebble => next
  end do
  if (not_done%n==0) then
    nullify(not_done%head, not_done%tail)
  end if
  if (verbose>3 .and. not_done%n>0) then
    print*,'gather: n, t_min, dt_min=', not_done%n, not_done%t_min, not_done%dt_min
  end if
  not_done%n_update = not_done%n_update + n
  call not_done%check_list ('after gather')
  if ((self%track .or. self%level>=track_level) .and. n>0) then
    print '(i5,1x,a,7i7,1p,2e12.3)', self%id, &
      ' gather: pebbles, not_done, accrete, remove, split, import, export, n', &
      self%pebbles%n, not_done%n_accrete, not_done%n_remove, not_done%n_split, &
      not_done%n_import, not_done%n_export, n, not_done%z_min
    self%pebbles%n_accrete = 0
    self%pebbles%n_export = 0
    self%pebbles%n_import = 0
    self%pebbles%n_remove = 0
    self%pebbles%n_frozen = 0
    self%pebbles%n_split = 0
  end if
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE gather

!===============================================================================
!> Bring all (remaining) pebbles up-to-date with respect to current time, by
!> running through update update sweeps until all pebbles are up-to-date.  We
!> first scatter all pebbles to a "not_done" buffer, and then gather them back
!> to the self%pebbles as they finish (some may be exported to neighbor ranks
!> in the process).
!===============================================================================
SUBROUTINE update_pebbles (self)
  class(pebble_patch_t):: self
  !.............................................................................
  class(nbor_t), pointer:: nbors, nbor
  class(*), pointer:: link
  class(pebble_t), pointer:: p, rm
  logical, save:: first_time=.true.
  integer:: n, verbose
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Add pebbles initially, and on the top and bottom surfaces
  !-----------------------------------------------------------------------------
  if (.not.self%do_pebbles) return
  call trace%begin ('pebbles_t%update_pebbles', itimer=itimer)
  if (self%pebbles%n>0) then
    if (verbose>1) then
      if (associated(self%pebbles%head)) then
        print'(a,2i8,1p,2e15.6)', 'update_pebbles: id, n, t, head%t=', self%id, self%pebbles%n, &
          self%time, self%pebbles%head%t
      else
        print *,'ERROR, no head, patch id', self%id,'  n should be', self%pebbles%n>0
        call trace%end (itimer)
        return
      end if
    end if
  else
    call trace%end (itimer)
    return
  end if
  !-----------------------------------------------------------------------------
  ! By using an anonymous task%link, which points back to the link that "owns"
  ! the task, we can get access the neighbor list.
  !-----------------------------------------------------------------------------
  link => self%link
  select type (link)
  type is (link_t)
  nbors => link%nbor
  end select
  !-----------------------------------------------------------------------------
  ! Now scatter the existing pebbles onto a "not_done" list, and gather them
  ! back as they reach the desired time.
  !-----------------------------------------------------------------------------
  if (verbose>2) &
    print*,'pebbles  in =', self%pebbles%n
  call self%scatter (nbors, self%not_done)
  if (verbose>3) &
    print*, 'not_done =', self%not_done%n, self%pebbles%dt_min, self%pebbles%dt_max
  n = 0
  verbose = 0
  !print *,self%id,'gather loop, n =', self%not_done%n, associated(self%not_done%head), &
    !self%not_done%n_remove, self%not_done%n_export
  do while (self%not_done%n > 0)
    if (verbose>3) &
      print*, 'not_done =', self%not_done%n, self%pebbles%dt_min, self%pebbles%dt_max
    call self%gather (self%not_done, nbors, verbose)
    n = n+1
    if (n > 10000) then
      print*, self%id, 'ERROR: gather loop stuck, removing pebbles:', self%not_done%n, associated(self%not_done%head)
      p => self%not_done%head
      do while (associated(p))
        print'(i8,4(2x,1p,3e10.2))', p%id, p%p, p%v, p%vg, p%t, p%dt, self%time
        rm => p
        p => p%next
        call self%not_done%remove (rm)
      end do
      self%not_done%n = 0
      verbose = 0
    else if (n > 9980) then
      verbose = 1
    end if
    if (debug) &
      print *,self%id,'gather loop, n =', self%not_done%n, associated(self%not_done%head), &
      self%not_done%n_remove, self%not_done%n_export
  end do
  if (verbose>3) &
    print*,'pebbles out =', self%pebbles%n, self%pebbles%dt_min, self%pebbles%dt_max
  call trace%end (itimer)
END SUBROUTINE update_pebbles

!===============================================================================
!> Print info to stdout
!===============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(pebble_patch_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !.............................................................................
  integer:: iv, nv
  !-----------------------------------------------------------------------------
  call self%gpatch_t%info (nq, ntask, experiment_name)
  if (io%verbose >= 0) then
    if (self%time>=self%print_next) then
      if (io%print_time > 0) &
        self%print_next = (nint(self%time/io%print_time)+1)*io%print_time
      nv = 1
      if (io%verbose>2) nv = self%nv
      if (io%verbose>0 .or. self%track) then
        if (self%nw > 1) then
          do iv=1,nv
            print '(a,i8,2i4,i3,f12.5,1p,2e12.4,2x,3(2e10.2,2x),z10,2i6)', &
              'upd', self%id, self%level, iv, omp%thread, self%time-self%dtime, &
              self%dtime, self%u_max, &
              minval(self%mem(:,:,:,iv,self%it ,2)), &
              maxval(self%mem(:,:,:,iv,self%it ,2)), &
              minval(self%mem(:,:,:,iv,self%it ,1)), &
              maxval(self%mem(:,:,:,iv,self%it ,1)), &
              minval(self%mem(:,:,:,iv,self%new,1)), &
              maxval(self%mem(:,:,:,iv,self%new,1)), &
              self%status, nq, ntask
          end do
        else
          do iv=1,nv
            print '(a,i8,2i4,i3,f12.5,1p,2e12.4,2x,2(2e14.6,2x),z10,2i6)', &
              'upd', self%id, self%level, iv, omp%thread, self%time-self%dtime, &
              self%dtime, self%u_max, &
              minval(self%mem(:,:,:,iv,self%it ,1)), &
              maxval(self%mem(:,:,:,iv,self%it ,1)), &
              minval(self%mem(:,:,:,iv,self%new,1)), &
              maxval(self%mem(:,:,:,iv,self%new,1)), &
              self%status, nq, ntask
          end do
        end if
      end if
    end if
  end if
  !if ((io%verbose>0.or.self%track) .and. self%pebbles%n>0) then
  !  print '(5(1x,i8,a,2x),1p,e11.2,2x,"t =",g12.3)', &
  !    self%id,               ' : has', &
  !    self%pebbles%n,        ' pebbles, exported', &
  !    self%pebbles%n_export, ', removed', &
  !    self%pebbles%n_remove, ', accreted', &
  !    self%pebbles%n_accrete,', min z', &
  !    self%pebbles%z_min, self%time
  !  self%pebbles%n_export = 0
  !end if
END SUBROUTINE info

!===============================================================================
!> Velocity interpolation. By using a 4-D pointer we can treat the velocity
!> component as a 3D+time array, which can be passed to interpolate_r4.
!===============================================================================
SUBROUTINE interpolate_vg (self, pebble, warn)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: pebble
  integer:: warn
  !.............................................................................
  real:: w(3), pos(3), v(3), pt(2)
  integer:: ii(3), j, k, jt(2)
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: velocity
  real, save:: h(3,3)=-1.0
  class(mesh_t), pointer:: m
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Cache staggered mesh index shift for the three velocity components
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('pebble_patch_t%interpolate_vg', itimer=itimer)
  if (h(1,1)==-1.0) then
    do j=1,3
      do k=1,3
        m => self%mesh(k)
        h(k,j) = m%h(self%idx%px+j-1)
      end do
    end do
    if (io%master) &
      print *,'interpolate_vg: h =',h
  end if
  !-----------------------------------------------------------------------------
  ! Add planet position and subtract patch position for interpolation
  !-----------------------------------------------------------------------------
  pos = pebble%p + (planet%position - self%position)
  warn = 0
  do j=1,3      ! velocity component
    do k=1,3    ! particle space coordinate
      m => self%mesh(k)
      w(k) = pos(k)/m%d - h(k,j) + m%o
      ii(k) = w(k)
      warn = warn + merge (1, 0, ii(k)<1 .or. ii(k)>m%gn-1)
      ii(k) = max(1,min(ii(k),m%gn-1))
      w(k) = w(k)-ii(k)
    end do
    velocity => self%vgas(:,:,:,j,:)
    pebble%vg(j) = self%interpolate_r4 (pebble%t, velocity, ii, w)
  end do
  if (pebble%id==id_dbg) then
    call time_interval (self, pebble%t, jt, pt)
    print '(a,1p,2e15.6,1x,2i3,1x,0p,2f5.2,1x,9i3)', 'time_interval: time,p%t,jt,pt,iit =', self%time, pebble%t, jt, pt, self%iit
    print '(a,1p,9e15.6)', 'time_interval: self%t(iit) =', self%t(self%iit)
    print '(a,i6,1x,3i4,2x,3f5.2,3i3,1x,2f5.2,1x,3f8.2,2(2x,3f7.3))', &
      'interpolate_vg:', self%id, ii, w, warn, jt, pt, pebble%p, pebble%v, pebble%vg
  end if
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE interpolate_vg

!===============================================================================
!===============================================================================
FUNCTION index_r4 (self, p, w, gn) RESULT (ii)
  class(pebble_patch_t):: self
  real, intent(in):: p(3)
  real, intent(out):: w(3)
  integer:: ii(3), gn(3)
  !-----------------------------------------------------------------------------
  w  = self%offset + (modulo(p-self%position+0.5_8*self%box,self%box)-0.5_8*self%box)/self%ds
  ii = w
  w  = w-ii
  if (any(ii<1 .or. ii>gn-1) .and. verbose>1) then
    !debug = .true.
    print*,self%id, 'index_r4: warning', ii, w, inside_r4 (self, p)
  end if
  ii = max(1,min(ii,gn-1))
  if (debug) then
    print *,'planet position:', planet%position
    print *,'  self position:', self%position
    print *,'pebble position:', p
    print *,'          index:', ii
    print *,'        weights:', w
  end if
END FUNCTION

!===============================================================================
!> Decide if a point is inside a patch, taking period wrapping into account
!===============================================================================
FUNCTION inside_r4 (self, p, debug)
  class(pebble_patch_t):: self
  real, intent(in):: p(3)
  logical:: inside_r4
  logical, optional:: debug
  real(8):: pos(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('inside_r4', itimer=itimer)
  pos = p - self%position
  where (self%periodic)
    pos = modulo(pos+self%box*0.5_8,self%box) - self%box*0.5_8
  end where
  inside_r4 = all(abs(pos) < 0.5*self%size)
  if (verbose>1 .and. present(debug)) then
    print *, self%id, 'remove:', p, real(self%box), real(abs(pos)-0.5*self%size)
  end if
  if (detailed_timer) call trace%end (itimer)
END FUNCTION

!=======================================================================
!=======================================================================
SUBROUTINE time_interval (self, time, jt, pt)
  class (pebble_patch_t)                 :: self
  real(8)                         :: time
  integer                         :: jt(2)
  real                            :: pt(2)
  !.....................................................................
  integer                         :: i, j
  real(8)                         :: t(self%nt), dt
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('pebble_patch_t%time_interval', 2, itimer=itimer)
  t = self%t(self%iit)
  j = 1
  do i=1,self%nt-2
    if (t(i) < time) j=i
  end do
  jt = self%iit(j:j+1)
  dt = t(j+1)-t(j)
  if (dt==0.0) then
    pt(2) = 0.0
  else
    pt(2) = (time-t(j))/dt
  end if
  pt(1) = 1.-pt(2)
  if (verbose>1 .and. pt(2)>0.0) then
    print '(a,1p,3e12.4,2i5,0p,2f8.4)', 'time_interval: times  =', t(j), time, t(j+1), jt, pt
  end if
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE time_interval

!===============================================================================
!===============================================================================
FUNCTION interpolate_r4 (self, t, field, i, w) RESULT (out)
  class(pebble_patch_t):: self
  real(8):: t
  real(kind=KindScalarVar), dimension(:,:,:,:), intent(in):: field
  integer:: i(3)
  real, intent(in):: w(3)
  real:: out
  !.............................................................................
  real:: q(3), pt(2)
  integer:: it, j, dim(4), jt(2)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call time_interval (self, t, jt, pt)
  if (detailed_timer) call trace%begin ('pebble_patch_t%interpolate_r4', itimer=itimer)
  q = 1.0-w
  dim = shape(field)-1
  if (any(i<1 .or. i>dim(1:3))) then
    if (verbose>1) print*,'interpolate_r4: warning', i, w
    i = max(1,min(i,dim(1:3)))
  end if
  out = &
    pt(1)*( &
      q(3)*(q(2)*(q(1)*field(i(1),i(2)  ,i(3)  ,jt(1))+w(1)*field(i(1)+1,i(2)  ,i(3)  ,jt(1)))  + &
            w(2)*(q(1)*field(i(1),i(2)+1,i(3)  ,jt(1))+w(1)*field(i(1)+1,i(2)+1,i(3)  ,jt(1)))) + &
      w(3)*(q(2)*(q(1)*field(i(1),i(2)  ,i(3)+1,jt(1))+w(1)*field(i(1)+1,i(2)  ,i(3)+1,jt(1)))  + &
            w(2)*(q(1)*field(i(1),i(2)+1,i(3)+1,jt(1))+w(1)*field(i(1)+1,i(2)+1,i(3)+1,jt(1))))) + &
    pt(2)*( &
      q(3)*(q(2)*(q(1)*field(i(1),i(2)  ,i(3)  ,jt(2))+w(1)*field(i(1)+1,i(2)  ,i(3)  ,jt(2)))  + &
            w(2)*(q(1)*field(i(1),i(2)+1,i(3)  ,jt(2))+w(1)*field(i(1)+1,i(2)+1,i(3)  ,jt(2)))) + &
      w(3)*(q(2)*(q(1)*field(i(1),i(2)  ,i(3)+1,jt(2))+w(1)*field(i(1)+1,i(2)  ,i(3)+1,jt(2)))  + &
            w(2)*(q(1)*field(i(1),i(2)+1,i(3)+1,jt(2))+w(1)*field(i(1)+1,i(2)+1,i(3)+1,jt(2)))))
  if (detailed_timer) call trace%end (itimer)
END FUNCTION

!===============================================================================
!===============================================================================
SUBROUTINE output (self)
  class(pebble_patch_t):: self
  !.............................................................................
  character(len=128):: filename
  integer, parameter:: data_format=5, n_var=9, n_id=3
  real, dimension(:,:), pointer:: buffer
  integer, dimension(:,:), pointer:: idbuf
  integer:: i, j, np, n_not
  class(pebble_t), pointer:: p
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.io%do_output .or. self%time < self%out_next) return
  !
  call self%gpatch_t%output
  call trace%begin ('pebble_patch_t%output', itimer=itimer)
  !
  if (self%do_pebbles) then
    !$omp critical (out_cr)
    write (filename,'(a,i5.5,"/",i5.5,".peb")') &
        trim(io%outputname), self%iout-1, self%id
    open  (io%data_unit, file=trim(filename), form='unformatted', status='unknown')
    write (io%data_unit) data_format
    !
    np = max(self%pebbles%n,1)
    allocate (buffer(np,n_var), idbuf(np,n_id))
    buffer = 0.0
    p => self%pebbles%head
    n_not = 0
    do i=1,np
      if (.not.associated(p)) then
        n_not = n_not+1
        cycle
      end if
      idbuf(i,1) = p%id
      idbuf(i,2) = p%status
      idbuf(i,3) = p%parent
      buffer(i,1:3) = p%p
      buffer(i,4:6) = p%v
      buffer(i,7:9) = [real(p%t,kind=4),p%s,p%w]
      p => p%next
    end do
    if (io%master .and. verbose>1) &
      print '(a,a,i7,a)', 'pebble_patch_t%output: ', trim(filename), np, ' particles'
    write (io%data_unit) np, n_var, shared%n_repeat, shared%n_patch, n_id
    write (io%data_unit) buffer
    write (io%data_unit) idbuf
    !
    close (io%data_unit)
    deallocate (buffer,idbuf)
    !$omp end critical (out_cr)
    if (self%id==2) call pebble_trace%flush_buffer                    ! Flush the pebble buffer every now and then
  end if
  call trace%end (itimer)
END SUBROUTINE output

!===============================================================================
!> Read .peb file produced by the output procedure above
!===============================================================================
SUBROUTINE input (self, ok)
  class(pebble_patch_t):: self
  logical:: ok
  !.............................................................................
  character(len=128):: filename
  real, dimension(:,:), pointer:: buffer
  integer, dimension(:,:), pointer:: idbuf
  integer:: i, np, n_repeat, n_patch, data_format, n_var
  class(pebble_t), pointer:: p
  logical:: exists
  !-----------------------------------------------------------------------------
  ok = .false.
  if (.not.self%do_pebbles.or.self%iout==0) return
  call trace%begin ('pebble_patch_t%input')
  !
  write (filename,'(a,"/",i5.5,"/",i5.5,".peb")') &
      trim(io%inputdir), self%iout-1, self%id
  inquire  (file=filename, exist=exists)
  if (exists) then
    if (verbose>1) print *, trim(filename), ' found'
  else
    if (verbose>0) print *, trim(filename), ' not found'
    call trace%end()
    return
  end if
  !
  !$omp critical (out_cr)
  open (io%data_unit, file=filename, form='unformatted', status='unknown')
  read (io%data_unit) data_format
  read (io%data_unit) np, n_var, n_repeat, n_patch
  allocate (buffer(np,n_var), idbuf(np,2))
  read (io%data_unit) buffer
  read (io%data_unit) idbuf
  if (io%master .and. verbose>0 .or. self%track) &
    print*, 'pebble_patch_t%input: ', trim(filename), np, ' particles'
  !
  do i=1,np
    allocate (p)
    p%id = idbuf(i,1)
    p%status = idbuf(i,2)
    p%p = buffer(i,1:3)
    p%v = buffer(i,4:6)
    p%t = buffer(i,7)
    p%s = buffer(i,8)
    p%w = buffer(i,9)
    p%m = p%s**3*self%pebbles%density
    call self%pebbles%append (p)
  end do
  close (io%data_unit)
  deallocate (buffer,idbuf)
  !$omp end critical (out_cr)
  ok = .true.
  call trace%end()
END SUBROUTINE input

END MODULE pebble_patch_mod
