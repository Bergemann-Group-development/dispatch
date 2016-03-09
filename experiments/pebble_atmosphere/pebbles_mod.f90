!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> A patch extension with pebbles.  There are so many patch properties that the
!> pebbles need, that it is much simpler to use an extension of an mhd_t, rather
!> than to pass (or/and recompute) those properties.  The key flag is do_pebbles,
!> which modifies slightly the variable allocations in mhd_mod.f90, so that
!> gas pressure and gas velocities are permanently allocated, and are available
!> to the pebble procedures.  To simplify neighbor access, the patch also has
!> an anonymous pointer to its primary list node, which has the nbor head of the
!> neighbor list.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pebbles_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE task_mod
  USE patch_mod
  USE mhd_mod
  USE list_mod
  USE planet_mod
  USE star_mod
  USE scaling_mod
  USE units_mod
  USE random_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! A single pebble in a linked list
  !-----------------------------------------------------------------------------
  type, public:: pebble_t
    integer:: id                        ! pebble nr
    integer:: it=0                      ! time step
    real::  p(3)                        ! position
    real::  v(3)                        ! velocity
    real:: vg(3)                        ! gas velocity
    real:: ag(3)                        ! acceleration of gravity
    real:: ad(3)                        ! drag acceleration
    real:: m                            ! mass
    real:: s                            ! size
    real:: d                            ! solid density
    real:: dg                           ! gas density
    real:: c                            ! drag acceleration rate
    real:: dt                           ! time step
    real(8):: t                         ! particle time
    type(pebble_t), pointer:: prev=>null(), next=>null()
  end type
  !-----------------------------------------------------------------------------
  ! A linked list of pebbles
  !-----------------------------------------------------------------------------
  type, public:: pebbles_t
    integer:: n                                 ! number of active pebbles
    integer:: n_old = 0                         ! have stepped to many times
    integer:: n_accrete = 0                     ! have accreted
    integer:: n_outside(3) = 0                  ! have left the patch
    integer:: n_behind = 0                      ! are behind self%time
    integer:: n_export = 0                      ! number of exported pebbles
    integer:: n_remove = 0                      ! number of removed pebbles
    integer:: it_max                            ! max life time
    integer:: id_dbg = 0                        ! follow this pebble
    integer:: ppc = 5                           ! particles per cell
    integer(8):: n_update = 0                   ! speed counter
    real:: mass                                 ! representative pebble mass
    real:: size                                 ! representative pebble size
    real:: density                              ! mass density
    real:: drag_coeff                           ! drag coefficient
    real:: t_min                                ! oldest t
    real:: dt_min, dt_max                       ! dt extrema
    real:: z_min                                ! smallest positive z
    real:: gamma                                ! gamma for sound speed
    real:: courant=0.1                          ! courant number for motion
    real:: position(3)                          ! relative to star
    type(pebble_t), pointer:: head, tail        ! first & last pebble
    logical:: initialized=.false.               ! initialize flag
  contains
    procedure append                            ! append a pebble
    procedure remove                            ! remove a pebble
    procedure export                            ! export to other patches (or remove)
    procedure delete                            ! delete pebbles
    procedure check_list
  end type
  !-----------------------------------------------------------------------------
  ! A patch with pebbles
  !-----------------------------------------------------------------------------
  type, public, extends(mhd_t):: pebble_patch_t
    class(pebbles_t), pointer:: pebbles=>null()
    class(pebbles_t), pointer:: not_done=>null()
    procedure(kdk_drag), pointer:: move=>null() ! particle mover
    real:: pebble_rate = 1e3                    ! pebble arrival rate
    real:: vz_pebble
  contains
    procedure init                              ! initialize the position and size
    procedure initial_add                       ! add pebbles
    procedure add                               ! add pebbles
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
  end type
  integer:: id=0
CONTAINS

!> =============================================================================
!> Create a number of particles proportional at the left hand side wall
!> =============================================================================
SUBROUTINE init (self)
  class(pebble_patch_t):: self
  !.............................................................................
  character(len=8), save:: type='kdk_drag'
  real, save:: density_yr_cgs = 1e-12, solid_density_cgs=3., size(3), pebble_size_cm=1.
  real, save:: pebble_rate_cgs = 1e3  ! about one pebble per second (3e7 per year)
  real, save:: vz_pebble_cgs = 3e6    ! about 30 km/s vertically
  integer, save:: it_max=1000000, id_dbg=0, ppc=5
  logical, save:: do_pebbles
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /pebble_params/ do_pebbles, type, density_yr_cgs, pebble_size_cm, &
    solid_density_cgs, size, it_max, pebble_rate_cgs, vz_pebble_cgs, id_dbg, &
    ppc
  !-----------------------------------------------------------------------------
  ! Make sure to initialize, but only once
  !-----------------------------------------------------------------------------
  if (.not.associated(self%pebbles)) allocate (self%pebbles)
  if (self%pebbles%initialized) return
  self%pebbles%initialized = .true.
  call trace_begin ('pebbles_t%init', 3)
  call random%init
  random%seed = 11
  call self%patch_t%init
  associate (mx=>self%mesh(1)%gn,my=>self%mesh(2)%gn,mz=>self%mesh(3)%gn,mt=>self%nt)
    allocate (self%density(mx,my,mz,mt))
    allocate (self%pressure(mx,my,mz,mt))
    allocate (self%vgas(mx,my,mz,3,mt))
  end associate
  !-----------------------------------------------------------------------------
  ! Read parameters, but only from one thread
  !-----------------------------------------------------------------------------
  !$omp barrier
  !$omp single
  if (first_time) then
    do_pebbles = .true.
    first_time = .false.
    size = self%size
    if (io%verbose>1) &
      print*, 'thread', omp_mythread(), ' reading namelist'
    rewind (io%input)
    read (io%input, pebble_params, iostat=iostat)
    if (io%master) write(*, pebble_params)
  end if
  !$omp end single
  !$omp barrier
  self%do_pebbles = do_pebbles
  self%pebble_rate = pebble_rate_cgs*scaling%t
  self%vz_pebble = vz_pebble_cgs/scaling%u
  !-----------------------------------------------------------------------------
  ! Compute drag coefficient; the inverse stopping time for a 1 cm particle in
  ! a gas density of 1 g/cm**3, in the sub-sonic Epstein regime, specified in
  ! terms of the cgs gas density that gives one year stopping time.
  !-----------------------------------------------------------------------------
  self%pebbles%ppc = ppc
  self%pebbles%id_dbg = id_dbg
  self%pebbles%size = pebble_size_cm/scaling%l
  self%pebbles%it_max = it_max
  self%pebbles%density = solid_density_cgs/scaling%d
  self%pebbles%mass = self%pebbles%density*self%pebbles%size**3
  self%pebbles%drag_coeff = scaling%d/(scaling%l*density_yr_cgs*cgs%yr/scaling%t)
  !-----------------------------------------------------------------------------
  ! Select mover
  !-----------------------------------------------------------------------------
  select case (type)
  case ('kdk_drag')
    self%move => kdk_drag
  end select
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Create a number of particles proportional at the left hand side wall,
!> with initially a velocity interpolated from the argument velocity field
!> =============================================================================
SUBROUTINE initial_add (self, n)
  class(pebble_patch_t):: self
  integer, intent(in), optional:: n
  !.............................................................................
  class(pebbles_t), pointer:: pebbles
  class(pebble_t), pointer:: pebble
  integer:: ix, iy, iz, n_local, iter, warn
  real:: w(3), sign
  integer, save:: nprint=10
  !-----------------------------------------------------------------------------
  if (.not.self%do_pebbles) return
  call trace_begin ('pebbles_t%add', 3)
  pebbles => self%pebbles
  if (present(n)) then
    n_local = n
  else
    n_local = self%pebble_rate*self%dtime
  end if
  if (io%verbose>1 .or. (io%verbose>0.and.nprint>0)) then
    nprint = nprint-1
    print *, ' ADDING PEBBLES TO PATCH', self%id
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
  do iz=self%mesh(3)%li,self%mesh(3)%ui
  do iy=self%mesh(2)%li,self%mesh(2)%ui
  do ix=self%mesh(1)%li,self%mesh(1)%ui
    do iter=1,self%pebbles%ppc
      !$omp atomic
      id = id+1
      allocate (pebble)
      pebble%id= id
      pebble%t = self%time
      pebble%s = self%pebbles%size*random%ran1()
      pebble%d = self%pebbles%density
      pebble%m = pebble%s**3*pebble%d
      pebble%p(1) =   (self%position(1) - planet%position(1)) + self%size(1)*(0.5-random%ran1())
      pebble%p(2) =   (self%position(2) - planet%position(2)) + self%size(2)*(0.5-random%ran1())
      pebble%p(3) =   (self%position(3) - planet%position(3)) + self%size(3)*(0.5-random%ran1())
      call interpolate_vg (self, pebble, warn)
      pebble%v = pebble%vg
      pebble%d = self%mem(ix,iy,iz,1,1,1)
      call self%pebbles%append (pebble)
    end do
  end do
  end do
  end do
  print *, 'initial_add:', self%id, self%position(3), self%pebbles%n
  !call self%pebbles%check_list ('after initial_add')
  call trace_end
END SUBROUTINE initial_add

!> =============================================================================
!> Create a number of particles proportional at the left hand side wall,
!> with initially a velocity interpolated from the argument velocity field
!> =============================================================================
SUBROUTINE add (self, n)
  class(pebble_patch_t):: self
  integer, intent(in), optional:: n
  !.............................................................................
  class(pebbles_t), pointer:: pebbles
  class(pebble_t), pointer:: pebble
  integer:: i, j, k, ii(3), n_local
  real:: w(3), sign
  integer, save:: nprint=10
  !-----------------------------------------------------------------------------
  if (.not.self%do_pebbles) return
  call trace_begin ('pebbles_t%add', 3)
  pebbles => self%pebbles
  if (present(n)) then
    n_local = n
  else
    n_local = self%pebble_rate*self%dtime
  end if
  if (io%verbose>1 .or. (io%verbose>0.and.nprint>0)) then
    nprint = nprint-1
    print *, ' ADDING PEBBLES TO PATCH', self%id
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
    !$omp atomic
    id = id+1
    allocate (pebble)
    pebble%id= id
    pebble%t = self%time
    pebble%s = self%pebbles%size
    pebble%d = self%pebbles%density
    pebble%m = pebble%s**3*pebble%d
    sign = merge(1.0,-1.0,self%position(3)>0.0)
    pebble%p(1) =   (self%position(1) - planet%position(1)) + self%size(1)*(0.5-random%ran1())
    pebble%p(2) =   (self%position(2) - planet%position(2)) + self%size(2)*(0.5-random%ran1())
    pebble%p(3) =   (self%position(3) - planet%position(3)) + self%size(3)*0.5*sign
    pebble%v(1) =   self%omega(3)*planet%orbit_radius
    pebble%v(2) =   0.0
    pebble%v(3) = - self%vz_pebble*sign
    call self%pebbles%append (pebble)
    self%pebbles%z_min = min(self%pebbles%z_min,abs(pebble%p(3)))
  end do
  !call self%pebbles%check_list ('after add')
  call trace_end
END SUBROUTINE add

!> =============================================================================
!> Create a number of particles proportional at the left hand side wall,
!> with initially a velocity interpolated from the argument velocity field
!> =============================================================================
SUBROUTINE delete (self)
  class(pebbles_t):: self
  !.............................................................................
  class(pebble_t), pointer:: pebble, old
  !-----------------------------------------------------------------------------
  call trace_begin ('pebbles_t%delete', 2)
  pebble => self%head
  do while (associated(pebble))
    old => pebble
    pebble => pebble%next
    deallocate (old)
  end do
  self%n = 0
  call trace_end
END SUBROUTINE delete

!> =============================================================================
!> =============================================================================
SUBROUTINE interpolate_vg (self, pebble, warn)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: pebble
  integer:: warn
  !.............................................................................
  real:: w(3), pos(3), v(3)
  integer:: ii(3), j, k
  real, dimension(:,:,:,:), pointer:: velocity
  !-----------------------------------------------------------------------------
  ! By using a 4-D pointer we can treat the velocity component as a 3D+time
  ! array, which can be passed to interpolate_r4.
  !-----------------------------------------------------------------------------
  pos = pebble%p + real(planet%position-self%position,kind=4)
  warn = 0
  do j=1,3
    do k=1,3
      w(k) = pos(k)/self%mesh(k)%d - self%mesh(k)%h(k+2) + self%mesh(k)%o
      ii(k) = w(k)
      warn = warn + merge (1, 0, ii(k)<1 .or. ii(k)>self%mesh(k)%gn-1)
      ii(k) = max(1,min(ii(k),self%mesh(k)%gn-1))
    end do
    velocity => self%vgas(:,:,:,j,:)
    pebble%vg(j) = self%interpolate_r4 (velocity, ii, w)
  end do
END SUBROUTINE interpolate_vg

!> =============================================================================
!> Look up gas velocity components, given positions, and compute the drag rate
!> (inverse stopping time)
!> =============================================================================
SUBROUTINE drag (self, pebble, warning)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: pebble
  !.............................................................................
  real:: w(3), pos(3), d, p, temp, gamma, v(3), speed2, sound2
  integer:: ii(3), j, k
  integer, save:: nprint=10
  integer:: w1, w2, warning
  !-----------------------------------------------------------------------------
  ii  = self%index_r4 (pebble%p, w, self%mesh%gn)
  d   = self%interpolate_r4 (self%density , ii, w)
  p   = self%interpolate_r4 (self%pressure, ii, w)
  w1  = count(w<0.0 .or. w>1.0)
  call interpolate_vg (self, pebble, w2)
  warning = w1+w2
  !
  speed2    = sum((pebble%v-pebble%vg)**2)
  temp      = p/d
  gamma     = merge(1.6666667,1.4,temp>2000.)
  sound2    = gamma*p/d
  pebble%c  = self%pebbles%drag_coeff*d/pebble%s*sqrt((speed2+sound2)/sound2)
  pebble%ad = (pebble%vg-pebble%v)*pebble%c
  if (w1+w2>0 .and. io%verbose>3) then
    pos = pebble%p + real(planet%position-self%position,kind=4)
    print '("ww",2i1,1p,9(3e10.2,1x))', &
      w1,w2, &
      pos, &
      pos+2.*pebble%dt*pebble%v, &
      pebble%v*pebble%dt/self%ds, &
      pebble%ag, pebble%dt, pebble%c*pebble%dt
  end if
  if (io%verbose>0 .and. nprint>0) then
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
    print *,'       time step (s):', pebble%dt*scaling%t
    print *,'           stiffness:', pebble%dt*pebble%c
    print *,'     cells traversed:', pebble%dt*pebble%v/self%ds
  end if
END SUBROUTINE drag

!> =============================================================================
!> Move a pebbles forward one step, with drag
!> =============================================================================
SUBROUTINE kdk_drag (self, p)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: p
  !.............................................................................
  real:: rp(3), ar
  integer, save:: nprint=20
  integer:: warn
  !-----------------------------------------------------------------------------
  ! Get acceleration from planet and Coriolis force from star, ag(t)
  !---------------------------------------------------------------------------
  rp = - p%p
  ar = 1./sqrt(sum(rp**2))
  p%ag = scaling%grav*planet%mass*rp*ar**3
  rp = self%pebbles%position - p%p
  ar = 1./sqrt(sum(rp**2))
  p%ag(1) = p%ag(1) + 2.*self%omega(3)*p%v(2)
  p%ag(2) = p%ag(2) - 2.*self%omega(3)*p%v(1) + self%omega(3)**2*p%p(2)
  p%ag(3) = p%ag(3)                           - self%omega(3)**2*p%p(3)
  !p%ag = p%ag + scaling%grav*star%mass*rp*ar**3
  p%dt = self%pebbles%courant*sqrt(sum(p%p**2)/sqrt(sum(p%ag**2)))
  p%dt = min(p%dt,minval(self%ds/max(1e-20,abs(p%v))))
  !
  if (p%id==self%pebbles%id_dbg) then
    print'("kdkd1",1p,4(3e11.2,2x))', p%p, p%v, p%ag, p%t
  end if
  !
  if (io%verbose>0 .and. p%id==1 .and. nprint>0) then
    nprint = nprint-1
    print'("kdkd:p,v,ag,ds,t,dt",1p,6(3e10.2,2x))', &
      p%p, p%v, p%ag, self%ds, p%t, p%dt
  end if
  !---------------------------------------------------------------------------
  ! Get the drag acceleration, at t+dt, then back up position
  !---------------------------------------------------------------------------
  p%p = p%p + p%dt*p%v
  p%t = p%t + p%dt
  call self%drag (p,warn)
  if (warn>0 .and. io%verbose>2) then
    print *,'  dt',p%dt,1
    print *,' dt1',self%pebbles%courant*sqrt(sum(p%v**2)/sum(p%ag**2))
    print *,' dt2',minval(self%ds/max(1e-20,abs(p%v)))
    print *,' dt3',min(p%dt,minval(self%ds/max(1e-20,abs(p%v))))
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
  ar = 1./sqrt(sum(rp**2))
  p%ag = scaling%grav*planet%mass*rp*ar**3
  rp = self%pebbles%position - p%p
  ar = 1./sqrt(sum(rp**2))
  !p%ag = p%ag + scaling%grav*star%mass*rp*ar**3
  p%ag(1) = p%ag(1) + 2.*self%omega(3)*p%v(2)
  p%ag(2) = p%ag(2) - 2.*self%omega(3)*p%v(1) + self%omega(3)**2*p%p(2)
  p%ag(3) = p%ag(3)                           - self%omega(3)**2*p%p(3)
  !---------------------------------------------------------------------------
  ! Get the drag acceleration, at t+2*dt
  !---------------------------------------------------------------------------
  call self%drag (p,warn)
  if (warn>0 .and. io%verbose>2) then
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
  p%v = p%v + p%dt*(p%ag+p%ad)/(1.+p%c*p%dt)
  if (p%id==self%pebbles%id_dbg .or. (io%verbose>0 .and. p%id==1 .and. nprint>0)) then
    nprint = nprint-1
    print'("kdkd2",1p,4(3e11.2,2x))', p%p, p%v, p%ag, p%t
  end if
END SUBROUTINE kdk_drag

!> =============================================================================
!> Check list consistency
!> =============================================================================
SUBROUTINE check_list (self, label)
  class(pebbles_t):: self
  character(len=*), optional:: label
  class(pebble_t), pointer:: p
  !.............................................................................
  integer:: n
  n=0
  p => self%head
  do while (associated(p))
    p => p%next
    n=n+1
  end do
  if (n /= self%n) then
    if (present(label)) then
      print*, 'inconsistent pebble list ', label, n, self%n
    else
      print*, 'inconsistent pebble list ', n, self%n
    end if
  end if
END SUBROUTINE check_list

!> =============================================================================
!> Append a pebble to the list
!> =============================================================================
SUBROUTINE append (self, pebble)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble
  !.............................................................................
  !$omp critical (pebble_cr)
  pebble%prev => self%tail
  if (.not.associated(self%head)) self%head => pebble
  if (associated(self%tail)) self%tail%next => pebble
  nullify (pebble%next)
  self%tail => pebble
  self%n = self%n + 1
  !$omp end critical (pebble_cr)
END SUBROUTINE append

!> =============================================================================
!> Remove a pebble from the list
!> =============================================================================
SUBROUTINE remove (self, pebble)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble, p
  !.............................................................................
  !$omp critical (pebble_cr)
  if (associated(pebble%prev)) then
    pebble%prev%next => pebble%next
  else
    self%head => pebble%next
    if (associated(self%head)) nullify(self%head%prev)
  end if
  if (associated(pebble%next)) then
    pebble%next%prev => pebble%prev
  else
    self%tail => pebble%prev
    if (associated(self%tail)) nullify(self%tail%next)
  end if
  self%n = self%n - 1
  !$omp end critical (pebble_cr)
END SUBROUTINE remove

!> =============================================================================
SUBROUTINE export (self, pebble, nbor_list)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble
  class(pebbles_t), pointer:: nbor_list
  integer, save:: nprint=20
  !.............................................................................
  call self%remove (pebble)
  if (associated(nbor_list)) then
    call nbor_list%append (pebble)
    self%n_export = self%n_export+1
  else
    if (nprint > 0) then
      nprint = nprint-1
      print *,'exporting pebble to a black hole ...'
    end if
  end if
  !call self%check_list ('after export')
END SUBROUTINE export

!> =============================================================================
!> =============================================================================
SUBROUTINE scatter (self, nbors, not_done)
  class(pebble_patch_t):: self
  class(link_t), pointer:: nbors
  !.............................................................................
  class(pebbles_t), pointer:: not_done
  class(pebble_t), pointer:: pebble, next
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  integer:: n
  !-----------------------------------------------------------------------------
  self%pebbles%dt_min = 1e30
  self%pebbles%dt_max = 0.0
  self%pebbles%z_min = 1e30
  n = 0
  pebble => self%pebbles%head
  do while (associated(pebble))
    n = n+1
    next => pebble%next
    call self%move (pebble)
    self%pebbles%z_min = min(self%pebbles%z_min,abs(pebble%p(3)))
    self%pebbles%dt_min = min(pebble%dt,self%pebbles%dt_min)
    self%pebbles%dt_max = max(pebble%dt,self%pebbles%dt_max)
    if (.not. self%inside_r4 (pebble%p+2.*pebble%dt*pebble%v)) then
      !-------------------------------------------------------------------------
      ! Check if it should be exported to a neighbor
      !-------------------------------------------------------------------------
      if (associated(nbors)) then
        nbor => nbors
        do while (associated(nbor))
          task => nbor%task
          select type (task)
          class is (pebble_patch_t)
            if (task%inside_r4(pebble%p)) then
              call self%pebbles%export (pebble, task%pebbles)
              goto 9
            end if
          end select
          nbor => nbor%next
        end do
        !-----------------------------------------------------------------------
        ! If not inside, and not inside a neighbor, kill it
        !-----------------------------------------------------------------------
        call self%pebbles%remove (pebble)
        self%pebbles%n_remove = self%pebbles%n_remove+1
        goto 9
      else
        !-----------------------------------------------------------------------
        ! No neighbors???
        !-----------------------------------------------------------------------
        print *,'WARNING: pebble_patch without neighbors', self%id
        call self%pebbles%remove (pebble)
        if (abs(pebble%p(1)) > self%size(1)*0.5) &
          self%pebbles%n_outside(1) = self%pebbles%n_outside(1) + 1
        if (abs(pebble%p(2)) > self%size(2)*0.5) &
          self%pebbles%n_outside(2) = self%pebbles%n_outside(2) + 1
        if (abs(pebble%p(3)) > self%size(3)*0.5) &
          self%pebbles%n_outside(3) = self%pebbles%n_outside(3) + 1
        goto 9
      end if
    else if (planet%inside(pebble%p)) then
      !-------------------------------------------------------------------------
      ! The pebble hit the planet; accreted
      !-------------------------------------------------------------------------
      call self%pebbles%remove (pebble)
      self%pebbles%n_accrete = self%pebbles%n_accrete + 1
      goto 9
    else if (pebble%t < self%time) then
      !-------------------------------------------------------------------------
      ! Finally, check if the pebble is inside a neighbor with higher level
      !-------------------------------------------------------------------------
      if (associated(nbors)) then
        nbor => nbors
        do while (associated(nbor))
          task => nbor%task
          if (task%level==self%level+1) then
            select type (task)
            class is (pebble_patch_t)
              if (task%inside_r4(pebble%p)) then
                call self%pebbles%export (pebble, task%pebbles)
                goto 9
              end if
            end select
          end if
          nbor => nbor%next
        end do
        !-----------------------------------------------------------------------
        ! If it wasn't, export to the 'not_done' list
        !-----------------------------------------------------------------------
        call self%pebbles%export (pebble, not_done)
      end if
    end if
  9 pebble => next
  end do
  !call self%pebbles%check_list ('after scatter')
  if (self%pebbles%n>0) then
    print*,'scatter: n, self%pebbles%n, not_done=', self%pebbles%n, n, self%not_done%n
    self%pebbles%n_update = self%pebbles%n_update + self%pebbles%n
  end if
END SUBROUTINE scatter

!> =============================================================================
!> =============================================================================
SUBROUTINE gather (self, pebbles, nbors, done)
  class(pebble_patch_t):: self
  class(pebbles_t), pointer:: pebbles
  class(link_t), pointer:: nbors
  !.............................................................................
  class(link_t), pointer:: nbor
  class(pebbles_t), pointer:: done
  class(pebble_t), pointer:: pebble
  class(task_t), pointer:: task
  integer:: n
  !-----------------------------------------------------------------------------
  pebbles%dt_min = 1e30
  pebbles%dt_max = 0.0
  pebbles%t_min = 1e30
  n = 0
  pebble => pebbles%head
  do while (associated(pebble))
    n = n+1
    call self%move (pebble)
    pebbles%t_min  = min(pebble%t ,pebbles%t_min)
    pebbles%dt_min = min(pebble%dt,pebbles%dt_min)
    pebbles%dt_max = max(pebble%dt,pebbles%dt_max)
    if (.not. self%inside_r4 (pebble%p+3.*pebble%dt*pebble%v)) then
      if (associated(nbors)) then
        nbor => nbors
        do while (associated(nbor))
          task => nbor%task
          select type (task)
          class is (pebble_patch_t)
            if (task%inside_r4(pebble%p)) then
              call pebbles%export (pebble, task%pebbles)
              goto 9
            end if
          nbor => nbor%next
          end select
        end do
        call pebbles%remove (pebble)
      else
        call pebbles%remove (pebble)
        if (abs(pebble%p(1)) > self%size(1)*0.5) &
          pebbles%n_outside(1) = pebbles%n_outside(1) + 1
        if (abs(pebble%p(2)) > self%size(2)*0.5) &
          pebbles%n_outside(2) = pebbles%n_outside(2) + 1
        if (abs(pebble%p(3)) > self%size(3)*0.5) &
          pebbles%n_outside(3) = pebbles%n_outside(3) + 1
      end if
    else if (planet%inside(pebble%p)) then
      if (io%verbose>1) &
        print*,'removing accreted pebble', pebble%id
      call pebbles%remove (pebble)
      pebbles%n_accrete = pebbles%n_accrete + 1
    else if (pebble%t >= self%time) then
      call pebbles%export (pebble, done)
    end if
  9 pebble => pebble%next
  end do
  if (pebbles%n>0) then
    print*,'gather: n, n, t_min=', pebbles%n, pebbles%t_min, pebbles%dt_min
    pebbles%n_update = pebbles%n_update + pebbles%n
  end if
  !call self%pebbles%check_list ('after gather')
END SUBROUTINE gather

!> =============================================================================
!> Bring all (remaining) pebbles up-to-date with respect to current time, by
!> running through update update sweeps until all pebbles are up-to-date.  We
!> first scatter all pebbles to a "not_done" buffer, and then gather them back
!> to the self%pebbles as they finish (some may be exported to neighbor ranks
!> in the process).
!> =============================================================================
SUBROUTINE update_pebbles (self)
  class(pebble_patch_t):: self
  !.............................................................................
  class(link_t), pointer:: nbors, nbor
  class(*), pointer:: link
  logical, save:: first_time
  !-----------------------------------------------------------------------------
  ! Add pebbles initially, and on the top and bottom surfaces
  !-----------------------------------------------------------------------------
  if (.not.self%do_pebbles) return
  if (self%time == 0d0) call self%initial_add
  if (self%level==1 .and. abs(self%position(3))>self%ds(3)) call self%add
  if (self%pebbles%n>0) then
    if (associated(self%pebbles%head)) then
      print*,'update_pebbles: id, n=', self%id, self%pebbles%n, &
        self%time, self%pebbles%head%t
    else
      print *,'ERROR, no head, patch id', self%id,'  n should be', self%pebbles%n>0
      return
    end if
  else
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
  if (first_time) then
    first_time = .false.
    nbor => nbors
    do while (associated(nbor))
      print*,self%id,' nbor:',nbor%task%id
      nbor => nbor%next
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Now scatter the existing pebbles onto a "not_done" list, and gather them
  ! back as they reach the desired time.
  !-----------------------------------------------------------------------------
  allocate (self%not_done)
  if (io%verbose>1) &
    print*,'pebbles  in =', self%pebbles%n
  call self%scatter (nbors, self%not_done)
  if (io%verbose>2) &
    print*, 'not_done =', self%not_done%n, self%pebbles%dt_min, self%pebbles%dt_max
  do while (self%not_done%n > 0)
    if (io%verbose>2) &
      print*, 'not_done =', self%not_done%n, self%pebbles%dt_min, self%pebbles%dt_max
    call self%gather (self%not_done, nbors, self%pebbles)
  end do
  if (io%verbose>1) &
    print*,'pebbles out =', self%pebbles%n, self%pebbles%dt_min, self%pebbles%dt_max
  deallocate (self%not_done)
END SUBROUTINE update_pebbles

!> =============================================================================
!> Print info to stdout
!> =============================================================================
SUBROUTINE info (self, nq, experiment_name)
  class(pebble_patch_t):: self
  integer, optional:: nq
  character(len=64), optional:: experiment_name
  !.............................................................................
  integer:: iv, nv
  !-----------------------------------------------------------------------------
  if (io%verbose >= 0) then
    if (self%time>=self%print_next) then
      if (self%print_time > 0) &
        self%print_next = (nint(self%time/self%print_time)+1)*self%print_time
      nv = 1
      if (io%verbose>2) nv = self%nv
      do iv=1,nv
        print '(a,i8,2i4,i3,f10.6,1p,e12.4,0p,f7.2,2x,1p,2(2e10.2,2x),2e14.6,i6,2x,2i3,2x,9i3)', &
          'upd',self%id, self%level, iv, omp_mythread(), self%time, self%dtime, self%courant, &
          minval(self%mem(:,:,:,iv,self%it ,2)), &
          maxval(self%mem(:,:,:,iv,self%it ,2)), &
          minval(self%mem(:,:,:,iv,self%it ,1)), &
          maxval(self%mem(:,:,:,iv,self%it ,1)), &
          minval(self%mem(:,:,:,iv,self%new,1)), &
          maxval(self%mem(:,:,:,iv,self%new,1)), &
          self%status, nq
      end do
    end if
  end if
  if (io%verbose>0 .and. self%pebbles%n>0) then
    print '(5(1x,i8,a,2x),1p,e11.2,0p,f11.3)', &
      self%id,               ' now has', &
      self%pebbles%n,        ' pebbles, exported', &
      self%pebbles%n_export, ', removed', &
      self%pebbles%n_remove, ', accreted', &
      self%pebbles%n_accrete,', min z', &
      self%pebbles%z_min, self%time
    self%pebbles%n_export = 0
  end if
END SUBROUTINE info

!> =============================================================================
!> =============================================================================
FUNCTION index_r4 (self, p, w, gn) RESULT (ii)
  class(pebble_patch_t):: self
  real, intent(in):: p(3)
  real, intent(out):: w(3)
  integer:: ii(3), gn(3)
  !-----------------------------------------------------------------------------
  w  = self%offset + (p+real(planet%position-self%position,kind=4))/self%ds
  ii = w
  w  = w-ii
  if (any(ii<1 .or. ii>gn-1) .and. io%verbose>1) then
    print*,'index_r4: warning', ii, w
    ii = max(1,min(ii,gn-1))
  end if
END FUNCTION

!> =============================================================================
!> =============================================================================
FUNCTION inside_r4 (self, p) RESULT (out)
  class(pebble_patch_t):: self
  real, intent(in):: p(3)
  logical:: out
  !-----------------------------------------------------------------------------
  out = all(abs(p+real(planet%position-self%position,kind=4)) < 0.5*self%size)
END FUNCTION

!> =============================================================================
!> =============================================================================
FUNCTION interpolate_r4 (self, field, i, w) RESULT (out)
  class(pebble_patch_t):: self
  real, dimension(:,:,:,:), intent(in):: field
  integer:: i(3)
  real, intent(in):: w(3)
  real:: out
  !.............................................................................
  real:: q(3)
  integer:: it, j, dim(4)
  !-----------------------------------------------------------------------------
  dim = shape(field)-1
  it = self%it  ! FIXME!
  q = 1.0-w
  if (any(i<1 .or. i>dim(1:3))) then
    if (io%verbose>1) print*,'interpolate_r4: warning', i, w
    i = max(1,min(i,dim(1:3)))
  end if
  out = &
    w(3)*(w(2)*(w(1)*field(i(1),i(2)  ,i(3)  ,it)+q(1)*field(i(1)+1,i(2)  ,i(3)  ,it))  + &
          q(2)*(w(1)*field(i(1),i(2)+1,i(3)  ,it)+q(1)*field(i(1)+1,i(2)+1,i(3)  ,it))) + &
    q(3)*(w(2)*(w(1)*field(i(1),i(2)  ,i(3)+1,it)+q(1)*field(i(1)+1,i(2)  ,i(3)+1,it))  + &
          q(2)*(w(1)*field(i(1),i(2)+1,i(3)+1,it)+q(1)*field(i(1)+1,i(2)+1,i(3)+1,it)))
END FUNCTION

!> =============================================================================
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(pebble_patch_t):: self
  character(len=64), optional:: experiment_name
  !.............................................................................
  character(len=128):: filename
  integer, parameter:: data_format=2, n_var=9
  real, dimension(:,:), pointer:: buffer
  integer, dimension(:), pointer:: idbuf
  integer:: i, j, np, n_not
  class(pebble_t), pointer:: p
  !-----------------------------------------------------------------------------
  if (.not.self%do_pebbles) return
  if (.not.io%do_output .or. self%time < self%out_next) return
  !$omp critical (pebble_cr)
  !
  !
  if (present(experiment_name)) then
    io%outputname = experiment_name
  end if
  if (trim(io%inputname) /= 'input') then
    io%outputname = io%inputname
  end if
  write (filename,'(a,i3.3,"_",i4.4,".peb")') &
      'data/'//trim(io%outputname)//'_', self%id, self%iout
  print *,'PEBBLE OUTPUT: ', trim(filename), self%pebbles%n
  open  (io%data_unit, file=filename, form='unformatted', status='unknown')
  write (io%data_unit) data_format
  !
  np = max(self%pebbles%n,1)
  allocate (buffer(np,n_var), idbuf(np))
  buffer = 0.0
  p => self%pebbles%head
  n_not = 0
  do i=1,np
    if (.not.associated(p)) then
      n_not = n_not+1
      cycle
    end if
    if (i<20) print*,p%id,p%p,p%v
    idbuf(i) = p%id
    buffer(i,1:3) = p%p
    buffer(i,4:6) = p%v
    buffer(i,7:9) = [real(p%t,kind=4),p%s,p%c]
    p => p%next
  end do
  print *,'PEBBLE OUTPUT 1', np, n_not
  write (io%data_unit) np, n_var
  write (io%data_unit) idbuf
  write (io%data_unit) buffer
  !
  close (io%data_unit)
  deallocate (buffer,idbuf)
  !
  call self%patch_t%output (experiment_name)
  !$omp end critical (pebble_cr)
END SUBROUTINE output

END MODULE pebbles_mod
