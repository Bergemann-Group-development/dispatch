!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> A self-contained data type for pebbles, to which one needs to  pass the
!> density, velocity, and patch information in argument lists
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
    integer:: it_max                            ! max life time
    integer(8):: n_update = 0                   ! speed counter
    real:: mass                                 ! representative pebble mass
    real:: size                                 ! representative pebble size
    real:: density                              ! mass density
    real:: drag_coeff                           ! drag coefficient
    real:: dt_min, dt_max                       ! dt extrema
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
  end type
  !-----------------------------------------------------------------------------
  ! A patch with pebbles
  !-----------------------------------------------------------------------------
  type, public, extends(mhd_t):: pebble_patch_t
    class(pebbles_t), pointer:: pebbles=>null()
    class(pebbles_t), pointer:: not_done=>null()
    procedure(kdk), pointer:: move=>null()      ! particle mover
  contains
    procedure init                              ! initialize the position and size
    procedure add                               ! add pebbles
    procedure drag                              ! compute gas drag
    procedure kdk                               ! kick-drift-kick, no drag
    procedure kdk_drag                          ! kick-drift-kick with drag
    procedure step_pebbles                      ! single step on all pebbles
    procedure scatter                           ! scatter to not_done
    procedure gather                            ! gather from not_done
    procedure update_pebbles                    ! advance until self%time
    procedure test                              ! run a test
    procedure index_r4                          ! compute interpolation index
    procedure inside_r4                         ! is it inside?
    procedure interpolate_r4                    ! interpolation
  end type
  integer:: id=0
CONTAINS

!===============================================================================
!> Create a number of particles proportional at the left hand side wall
!===============================================================================
SUBROUTINE init (self)
  class(pebble_patch_t):: self
  !.............................................................................
  character(len=8), save:: type='kdk'
  real:: density_yr_cgs = 1e-12, solid_density_cgs=3., size(3), pebble_size_cm=1.
  integer:: it_max=1000000
  logical:: do_pebbles
  integer:: iostat
  namelist /pebble_params/ do_pebbles, type, density_yr_cgs, pebble_size_cm, &
    solid_density_cgs, size, it_max
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
  do_pebbles = .true.
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
  size = self%size
  if (io%verbose>1) &
    print*, 'thread', omp_mythread(), ' reading namelist'
  rewind (io%input)
  read (io%input, pebble_params, iostat=iostat)
  if (io%master) write(*, pebble_params)
  !$omp end single
  !$omp barrier
  self%do_pebbles = do_pebbles
  !-----------------------------------------------------------------------------
  ! Compute drag coefficient; the inverse stopping time for a 1 cm particle in
  ! a gas density of 1 g/cm**3, in the sub-sonic Epstein regime, specified in
  ! terms of the cgs gas density that gives one year stopping time.
  !-----------------------------------------------------------------------------
  self%pebbles%size = pebble_size_cm/scaling%l
  self%pebbles%it_max = it_max
  self%pebbles%density = solid_density_cgs/scaling%d
  self%pebbles%mass = self%pebbles%density*self%pebbles%size**3
  self%pebbles%drag_coeff = scaling%d/(scaling%l*density_yr_cgs*cgs%yr/scaling%t)
  !-----------------------------------------------------------------------------
  ! Select mover
  !-----------------------------------------------------------------------------
  select case (type)
  case ('kdk')
    self%move => kdk
  case ('kdk_drag')
    self%move => kdk_drag
  end select
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Create a number of particles proportional at the left hand side wall,
!> with initially a velocity interpolated from the argument velocity field
!===============================================================================
SUBROUTINE add (self, n)
  class(pebble_patch_t):: self
  integer, intent(in):: n
  !.............................................................................
  class(pebbles_t), pointer:: pebbles
  class(pebble_t), pointer:: pebble
  integer:: i, j, k, ii(3)
  real:: w(3)
  !-----------------------------------------------------------------------------
  call trace_begin ('pebbles_t%add', 3)
  pebbles => self%pebbles
  if (io%master) then
    print 1, '   pebble size (cm):'   , self%pebbles%size*scaling%l
    print 1, 'solid density (cgs):', self%pebbles%density*scaling%d
    print 1, '         mass (cgs):', self%pebbles%mass*scaling%m
    print 1, '   model size (cgs):', self%size*scaling%l
    print 1, ' gas velocity (cgs):', self%vgas(3,3,3,1:3,self%it)*scaling%u
    1 format(a,1p,3e11.2)
  end if
  do i=1,n
    !$omp atomic
    id = id+1
    allocate (pebble)
    pebble%id= id
    pebble%t = self%time
    pebble%s = self%pebbles%size
    pebble%d = self%pebbles%density
    pebble%m = pebble%s**3*pebble%d
    pebble%p(1) = - self%size(1)*0.5
    pebble%p(2) =   self%size(2)*(0.5-random%ran1())
    pebble%p(3) =   self%size(3)*(0.5-random%ran1())
    block
      real, dimension(:,:,:,:), pointer:: velocity
      do j=1,3
        do k=1,3
          w(k) = (pebble%p(k) - self%mesh(k)%h(k+2))/self%mesh(k)%d + self%mesh(k)%o
          ii(k) = w(k)
          w(k) = w(k) - ii(k)
        end do
        velocity => self%vgas(:,:,:,j,:)
        pebble%vg(j) = self%interpolate_r4 (velocity, ii, w)
      end do
    end block
    pebble%v  = pebble%vg
    if (io%verbose>2) &
      print'(i7,1p,e10.2,e10.2,2(2x,3e10.2))', pebble%id, pebble%s, pebble%m, pebble%p, pebble%v
    call self%pebbles%append (pebble)
  end do
  call trace_end
END SUBROUTINE add

!===============================================================================
!> Create a number of particles proportional at the left hand side wall,
!> with initially a velocity interpolated from the argument velocity field
!===============================================================================
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

!===============================================================================
!> Kick step of kick-drift-kick (KDK), with gravity from star and planet
!===============================================================================
SUBROUTINE kdk (self, p)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: p
  !.............................................................................
  real:: ar, rp(3)
  integer, save:: nprint=20
  !-----------------------------------------------------------------------------
  rp = - p%p
  ar = 1./sqrt(sum(rp**2))
  p%ag = scaling%grav*planet%mass*rp*ar**3
  rp = self%pebbles%position - p%p
  ar = 1./sqrt(sum(rp**2))
  !p%ag = p%ag + scaling%grav*star%mass*rp*ar**3
  p%ag(1) = p%ag(1) + 2.*self%omega(3)*p%v(2)
  p%ag(2) = p%ag(2) - 2.*self%omega(3)*p%v(1) + self%omega(3)**2*p%p(2)
  p%dt = self%pebbles%courant*sqrt(sum(p%v**2)/sum(p%ag**2))
  p%dt = min(p%dt,minval(self%ds/max(1e-20,abs(p%v))))
  if (io%verbose>0 .and. p%id==1 .and. nprint>0) then
    nprint = nprint-1
    print'(" kdk",1p,4(3e11.2,2x))', p%p, p%v, p%ag, p%t, p%dt
  end if
  p%v = p%v + p%dt*p%ag

  p%p = p%p + 2.*p%dt*p%v
  p%t = p%t + 2.*p%dt

  rp = - p%p
  ar = 1./sqrt(sum(rp**2))
  p%ag = scaling%grav*planet%mass*rp*ar**3
  rp = self%pebbles%position - p%p
  ar = 1./sqrt(sum(rp**2))
  !p%ag = p%ag + scaling%grav*star%mass*rp*ar**3
  p%ag(1) = p%ag(1) + 2.*self%omega(3)*p%v(2)
  p%ag(2) = p%ag(2) - 2.*self%omega(3)*p%v(1) + self%omega(3)**2*p%p(2)
  p%v = p%v + p%dt*p%ag
  if (io%verbose>0 .and. p%id==1 .and. nprint>0) then
    nprint = nprint-1
    print'(" kdk",1p,4(3e11.2,2x))', p%p, p%v, p%ag, p%t
  end if
END SUBROUTINE kdk

!===============================================================================
!> Look up gas velocity components, given positions, and compute the drag rate
!> (inverse stopping time)
!===============================================================================
SUBROUTINE drag (self, pebble, warning)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: pebble
  !.............................................................................
  real:: w(3), d, p, temp, gamma, v(3), speed2, sound2
  integer:: ii(3), j, k
  integer, save:: nprint=10
  integer:: w1, w2, warning
  !-----------------------------------------------------------------------------
  ii        = self%index_r4 (pebble%p, w, self%mesh%gn)
  d         = self%interpolate_r4 (self%density , ii, w)
  p         = self%interpolate_r4 (self%pressure, ii, w)
  w1 = count(w<0.0 .or. w>1.0)
  !-----------------------------------------------------------------------------
  ! By using a 4-D pointer we can treat the velocity component as a 3D+time
  ! array, which can be passed to interpolate_r4.
  !-----------------------------------------------------------------------------
  block
    real, dimension(:,:,:,:), pointer:: velocity
    w2 = 0
    do j=1,3
      do k=1,3
        w(k) = pebble%p(k)/self%mesh(k)%d - self%mesh(k)%h(k+2) + self%mesh(k)%o
        ii(k) = w(k)
        ii(k) = max(1,min(ii(k),self%mesh(k)%gn))
        w(k) = w(k) - ii(k)
        if (w(k)<0.0 .or. w(k)>1.0) then
          print '("w:p,p1,p2,v*dt/ds,h,o=",1p,10e12.3))', &
            pebble%p(k), &
            pebble%p(k)-pebble%dt*pebble%v(k), &
            pebble%p(k)-pebble%dt*pebble%v(k)*2., &
            pebble%v(k)*pebble%dt/self%ds(k), &
            self%mesh(k)%h(k+2), self%mesh(k)%o
        end if
      end do
      w2 = w2 + count(w<0.0 .or. w>1.0)
    end do
  end block
  warning = w1+w2
  speed2    = sum((pebble%v-pebble%vg)**2)
  temp      = p/d
  gamma     = merge(1.6666667,1.4,temp>2000.)
  sound2    = gamma*p/d
  pebble%c  = self%pebbles%drag_coeff*d/pebble%s*sqrt((speed2+sound2)/sound2)
  pebble%ad = (pebble%vg-pebble%v)*pebble%c
  if (w1+w2>0) &
    print '("w",2i1,1p,9(3e10.2,1x))', &
      w1,w2, &
      pebble%p, &
      pebble%p+2.*pebble%dt*pebble%v, &
      pebble%v*pebble%dt/self%ds, &
      pebble%ag, pebble%dt, pebble%c*pebble%dt
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
  end if
END SUBROUTINE drag

!===============================================================================
!> Move a pebbles forward one step, with drag
!===============================================================================
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
  !p%ag = p%ag + scaling%grav*star%mass*rp*ar**3
  p%ag(1) = p%ag(1) + 2.*self%omega(3)*p%v(2)
  p%ag(2) = p%ag(2) - 2.*self%omega(3)*p%v(1) + self%omega(3)**2*p%p(2)
  p%dt = self%pebbles%courant*sqrt(sum(p%v**2)/sum(p%ag**2))
  p%dt = min(p%dt,minval(self%ds/max(1e-20,abs(p%v))))
  if (io%verbose>0 .and. p%id==1 .and. nprint>0) then
    nprint = nprint-1
    print'("kdkd",1p,4(3e11.2,2x))', p%p, p%v, p%ag, p%t, p%dt
  end if

  !---------------------------------------------------------------------------
  ! Get the drag acceleration, at t+dt, then back up position
  !---------------------------------------------------------------------------
  p%p = p%p + p%dt*p%v
  p%t = p%t + p%dt
  call self%drag (p,warn)
  if (warn>0) then
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

  !---------------------------------------------------------------------------
  ! Get the drag acceleration, at t+2*dt
  !---------------------------------------------------------------------------
  call self%drag (p,warn)
  if (warn>0) then
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
  if (io%verbose>0 .and. p%id==1 .and. nprint>0) then
    nprint = nprint-1
    print'("kdkd",1p,4(3e11.2,2x))', p%p, p%v, p%ag, p%t
  end if
END SUBROUTINE kdk_drag

!===============================================================================
!> Append a pebble to the list
!===============================================================================
SUBROUTINE append (self, pebble)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble
  !.............................................................................
  pebble%prev => self%tail
  if (.not.associated(self%head)) self%head => pebble
  if (associated(self%tail)) self%tail%next => pebble
  nullify (pebble%next)
  self%tail => pebble
  self%n = self%n + 1
END SUBROUTINE append

!===============================================================================
!> Remove a pebble from the list
!===============================================================================
SUBROUTINE remove (self, pebble)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble, p
  !.............................................................................
  if (associated(pebble%next)) then
    pebble%next%prev => pebble%prev
  else
    self%tail => pebble%prev
  end if
  if (associated(pebble%prev)) then
    pebble%prev%next => pebble%next
  else
    self%head => pebble%next
  end if
  self%n = self%n - 1
END SUBROUTINE remove

!===============================================================================
SUBROUTINE export (self, pebble, nbor_list)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble
  class(pebbles_t), pointer:: nbor_list
  integer, save:: nprint=20
  !.............................................................................
  call self%remove (pebble)
  if (associated(nbor_list)) then
    call nbor_list%append (pebble)
  else
    if (nprint > 0) then
      nprint = nprint-1
      print *,'exporting pebble to a black hole ...'
    end if
  end if
END SUBROUTINE export

!===============================================================================
!> Move the pebbles forward one sweep, and check if still inside our patch.
!> If not, search our neighbor list, and export to the neighbor pebble list.
!> Also check if it has hit the planet, and count that as an accretion
!===============================================================================
SUBROUTINE step_pebbles (self, nbors)
  class(pebble_patch_t):: self
  class(link_t), pointer:: nbors
  !.............................................................................
  class(link_t), pointer:: nbor
  class(pebble_t), pointer:: pebble
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  pebble => self%pebbles%head
  self%pebbles%dt_min = 1e30
  self%pebbles%dt_max = 0.0
  do while (associated(pebble))
    call self%move (pebble)
    self%pebbles%dt_min = min(pebble%dt,self%pebbles%dt_min)
    self%pebbles%dt_max = max(pebble%dt,self%pebbles%dt_max)
    if (pebble%it > self%pebbles%it_max) then
      call self%pebbles%remove (pebble)
      self%pebbles%n_old = self%pebbles%n_old + 1
    else if (.not. self%inside_r4 (pebble%p+pebble%dt*pebble%v)) then
      if (associated(nbors)) then
        nbor => nbors
        do while (associated(nbor))
          task => nbor%task
          print*,'NBOR!'
          stop
          select type (task)
          class is (pebble_patch_t)
            if (task%inside_r4(pebble%p)) then
              call self%pebbles%export (pebble, task%pebbles)
              exit
            end if
          class is (patch_t)
            call self%pebbles%remove (pebble)
            exit
          end select
          nbor => nbor%next
        end do
      else
        call self%pebbles%remove (pebble)
        if (abs(pebble%p(1)) > self%size(1)*0.5) &
          self%pebbles%n_outside(1) = self%pebbles%n_outside(1) + 1
        if (abs(pebble%p(2)) > self%size(2)*0.5) &
          self%pebbles%n_outside(2) = self%pebbles%n_outside(2) + 1
        if (abs(pebble%p(3)) > self%size(3)*0.5) &
          self%pebbles%n_outside(3) = self%pebbles%n_outside(3) + 1
      end if
    else if (planet%inside(pebble%p)) then
      call self%pebbles%remove (pebble)
      self%pebbles%n_accrete = self%pebbles%n_accrete + 1
    end if
    pebble%it = pebble%it+1
    if (pebble%t > self%time) then
      call self%pebbles%remove (pebble)
    end if
    pebble => pebble%next
  end do
  self%pebbles%n_update = self%pebbles%n_update + self%pebbles%n
END SUBROUTINE step_pebbles

!===============================================================================
!===============================================================================
SUBROUTINE scatter (self, pebbles, nbors, not_done)
  class(pebble_patch_t):: self
  class(pebbles_t), pointer:: pebbles
  class(link_t), pointer:: nbors
  !.............................................................................
  class(link_t), pointer:: nbor
  class(pebbles_t), pointer:: not_done
  class(pebble_t), pointer:: pebble
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  pebble => pebbles%head
  self%pebbles%dt_min = 1e30
  self%pebbles%dt_max = 0.0
  do while (associated(pebble))
    call self%move (pebble)
    self%pebbles%dt_min = min(pebble%dt,self%pebbles%dt_min)
    self%pebbles%dt_max = max(pebble%dt,self%pebbles%dt_max)
    if (.not. self%inside_r4 (pebble%p+2.*pebble%dt*pebble%v)) then
      if (associated(nbors)) then
        nbor => nbors
        do while (associated(nbor))
          task => nbor%task
          select type (task)
          class is (pebble_patch_t)
            if (task%inside_r4(pebble%p)) then
              call pebbles%export (pebble, task%pebbles)
              exit
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
      call pebbles%remove (pebble)
      pebbles%n_accrete = pebbles%n_accrete + 1
    else if (pebble%t < self%time) then
      call pebbles%export (pebble, not_done)
    end if
    pebble => pebble%next
  end do
  pebbles%n_update = pebbles%n_update + pebbles%n
END SUBROUTINE scatter

!===============================================================================
!===============================================================================
SUBROUTINE gather (self, pebbles, nbors, done)
  class(pebble_patch_t):: self
  class(pebbles_t), pointer:: pebbles
  class(link_t), pointer:: nbors
  !.............................................................................
  class(link_t), pointer:: nbor
  class(pebbles_t), pointer:: done
  class(pebble_t), pointer:: pebble
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  pebble => pebbles%head
  self%pebbles%dt_min = 1e30
  self%pebbles%dt_max = 0.0
  do while (associated(pebble))
    call self%move (pebble)
    self%pebbles%dt_min = min(pebble%dt,self%pebbles%dt_min)
    self%pebbles%dt_max = max(pebble%dt,self%pebbles%dt_max)
    if (.not. self%inside_r4 (pebble%p+3.*pebble%dt*pebble%v)) then
      if (associated(nbors)) then
        nbor => nbors
        do while (associated(nbor))
          task => nbor%task
          select type (task)
          class is (pebble_patch_t)
            if (task%inside_r4(pebble%p)) then
              call pebbles%export (pebble, task%pebbles)
              exit
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
    pebble => pebble%next
  end do
  pebbles%n_update = pebbles%n_update + pebbles%n
END SUBROUTINE gather

!===============================================================================
!> Bring all (remaining) pebbles up-to-date with respect to current time, by
!> running through update update sweeps until all pebbles are up-to-date.  We
!> first scatter all pebbles to a "not_done" buffer, and then gather them back
!> to the self%pebbles as they finish (some may be exported to neighbor ranks
!> in the process).
!===============================================================================
SUBROUTINE update_pebbles (self, nbors)
  class (pebble_patch_t):: self
  class (link_t), pointer:: nbors
  !.............................................................................
  allocate (self%not_done)
  if (io%verbose>1) &
    print*,'pebbles  in =', self%pebbles%n
  call self%scatter (self%pebbles, nbors, self%not_done)
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
  
!===============================================================================
!===============================================================================
FUNCTION index_r4 (self, p, w, gn) RESULT (ii)
  class(pebble_patch_t):: self
  real, intent(in):: p(3)
  real, intent(out):: w(3)
  integer:: ii(3), gn(3)
  !-----------------------------------------------------------------------------
  w  = self%offset + p/self%ds
  ii = w
  ii = max(1,min(ii,gn))
  w  = w-ii
END FUNCTION

!===============================================================================
!===============================================================================
FUNCTION inside_r4 (self, p) RESULT (out)
  class(pebble_patch_t):: self
  real, intent(in):: p(3)
  logical:: out
  !-----------------------------------------------------------------------------
  out = all(abs(p) < 0.5*self%size)
END FUNCTION

!===============================================================================
!===============================================================================
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
  dim = shape(field)
  it = self%it  ! FIXME!
  q = 1.0-w
  if (any(i<1 .or. i>dim(1:3))) then
    print*,'interpolate_r4: warning', i, w
    i = max([1,1,1],min(i,dim(1:3)))
  end if
  out = &
    w(3)*(w(2)*(w(1)*field(i(1),i(2)  ,i(3)  ,it)+q(1)*field(i(1)+1,i(2)  ,i(3)  ,it))  + &
          q(2)*(w(1)*field(i(1),i(2)+1,i(3)  ,it)+q(1)*field(i(1)+1,i(2)+1,i(3)  ,it))) + &
    q(3)*(w(2)*(w(1)*field(i(1),i(2)  ,i(3)+1,it)+q(1)*field(i(1)+1,i(2)  ,i(3)+1,it))  + &
          q(2)*(w(1)*field(i(1),i(2)+1,i(3)+1,it)+q(1)*field(i(1)+1,i(2)+1,i(3)+1,it)))
END FUNCTION

!===============================================================================
!> Create a number of particles and move them with an artificial gas density
!> and velocity field.  The calling sequence is
!> 
!> PROGRAM
!>   allocate (patch)           ! sets patch%size=1, patch%n=
!>   pebble_mod::test           !
!>     scaling::init
!>     pebble_mod::init
!>       patch_mod::init
!>         mesh_mod::init
!===============================================================================
SUBROUTINE test (self, link)
  class(pebble_patch_t):: self
  !.............................................................................
  class(pebbles_t), pointer:: pebbles
  class(pebble_t), pointer:: pebble
  class(link_t), pointer:: link
  real:: position(3), size(3)
  real:: gas_density_cgs=1e-10, gas_temperature_k=1000.
  integer:: k, it, ignore, it_max=2**20, niter, n_pebbles=10000
  namelist /test_params/ gas_density_cgs, gas_temperature_k, it_max, n_pebbles
  !-----------------------------------------------------------------------------
  call trace_begin ('pebbles_t%test')
  rewind (io%input)
  read (io%input, test_params, iostat=ignore)
  if (io%master) write (*, test_params)
  !-----------------------------------------------------------------------------
  ! Specify a star and a planet, and a 500 R_Earth box
  !-----------------------------------------------------------------------------
  call scaling%init
  star%mass   = cgs%m_sun/scaling%m
  planet%mass = cgs%m_earth/scaling%m
  star%position = 0.0
  planet%position = [0.0,1.0,0.0]*cgs%au/scaling%l
  planet%radius = cgs%r_earth/scaling%l
  self%omega(3) = 2.*cgs%pi/(cgs%yr/scaling%t)
  size = 500*planet%radius
  self%size = size
  self%position = planet%position
  !-----------------------------------------------------------------------------
  ! By calling the init hierarchy with do_pebbles true, we get access to the
  ! velocities and densities at all time slices, as selv%density and self%vgas
  !-----------------------------------------------------------------------------
  self%do_pebbles = .true.
  call self%init
  self%density => self%mem(:,:,:,1,:,1)
  self%pebbles%position = star%position - planet%position
  self%pebbles%it_max = it_max
  !-----------------------------------------------------------------------------
  ! Gas properties
  !-----------------------------------------------------------------------------
  self%density = gas_density_cgs/scaling%d
  self%pressure = self%density*gas_temperature_k/scaling%temp
  do it=1,self%nt
    self%vgas(:,:,:,1,it) = 0.5*cgs%kms/scaling%u*it
    self%vgas(:,:,:,2,it) = 0.0
    self%vgas(:,:,:,3,it) = 0.0
  end do
  if (io%verbose>1) then
    block
    real, dimension(:,:,:,:), pointer:: velocity
    do it=1,self%nt
      velocity => self%vgas(:,:,:,1,:)
      print*,velocity(3,3,3,it)
      print*,associated(self%vgas)
      print*, 'vg_min,max', minval(self%vgas(:,:,:,1,it)), maxval(self%vgas(:,:,:,1,it))
      print*, 've_min,max', minval(velocity(:,:,:,it)), maxval(velocity(:,:,:,it))
    end do
    end block
  end if
  if (io%master.and.io%verbose>0) then
    print *,'PEBBLE_PATCH PROPERTIES'
    print *,'position:', self%position
    print *,'position:', self%position
    print *,'    size:', self%size
    print *,'      ds:', self%ds
    print *,'       n:', self%n
    print *,'      ng:', self%ng
    print *,'      nv:', self%nv
    print *,'      nt:', self%nt
  end if
  !-----------------------------------------------------------------------------
  ! Test the simple KDK
  !-----------------------------------------------------------------------------
  random%seed = 11
  self%time = 0.0
  call self%add (n_pebbles)
  self%move => kdk
  self%dtime = self%ds(1)/2.
  niter = 50
  call tic
  do it=1,niter
    self%time = self%time + self%dtime
    call self%update_pebbles (link%nbor)
  end do
  call timer%toc ('kdk', self%pebbles%n_update)
  print 1, &
    'KDK   accreted:', self%pebbles%n_accrete, self%pebbles%n_outside, self%pebbles%n_old
  1 format(1x,a,i7,10x,"walls:",3i7,10x,"old:",i7)
  !-----------------------------------------------------------------------------
  ! Test the KDK + drag, with a 0.5 km/s head wind, blowing left to right
  !-----------------------------------------------------------------------------
  self%pebbles%n_update = 0
  self%pebbles%n_accrete = 0
  self%pebbles%n_outside = 0
  self%pebbles%n_old = 0
  random%seed = 11
  id = 0
  self%time = 0.0
  call self%pebbles%delete
  call self%add (n_pebbles)
  self%move => kdk_drag
  call tic
  do it=1,niter
    self%time = self%time + self%dtime
    call self%update_pebbles (link%nbor)
  end do
  call timer%toc ('kdk_drag', self%pebbles%n_update)
  print 1, &
    'KDK_d accreted:', self%pebbles%n_accrete, self%pebbles%n_outside, self%pebbles%n_old
  call trace_end
END SUBROUTINE test

END MODULE pebbles_mod
