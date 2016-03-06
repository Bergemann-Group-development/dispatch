!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> A self-contained data type for pebbles, to which one needs to  pass the
!> density, velocity, and patch information in argument lists
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pebbles_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE task_mod
  USE patch_mod
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
    type(pebble_t), pointer:: prev, next
  end type
  !-----------------------------------------------------------------------------
  ! A linked list of pebbles
  !-----------------------------------------------------------------------------
  type, public:: pebbles_t
    integer:: n                                 ! number of active pebbles
    integer:: n_accrete = 0                     ! has accreted
    integer(8):: n_update = 0                   ! speed counter
    real:: mass                                 ! representative pebble mass
    real:: size                                 ! representative pebble size
    real:: density                              ! mass density
    real:: drag_coeff                           ! drag coefficient
    real:: gamma                                ! gamma for sound speed
    real:: courant=0.1                          ! courant number for motion
    real:: t_max                                ! max life time
    real:: position(3)                          ! relative to star
    type(pebble_t), pointer:: head, tail        ! first & last pebble
    logical:: initialized=.false.               ! initialize flag
  contains
    procedure append                            ! append a pebble
    procedure remove                            ! remove a pebble
  end type
  !-----------------------------------------------------------------------------
  ! A patch with pebbles
  !-----------------------------------------------------------------------------
  type, public, extends(patch_t):: pebble_patch_t
    class(pebbles_t), pointer:: pebbles=>null()
    procedure(kdk), pointer:: move=>null()      ! particle mover
  contains
    procedure init                              ! initialize the position and size
    procedure add                               ! add pebbles
    procedure drag                              ! compute gas drag
    procedure kdk                               ! kick-drift-kick, no drag
    procedure kdk_drag                          ! kick-drift-kick with drag
    procedure update_pebbles
    procedure export                            ! export to other patches (or remove)
    procedure test                              ! run a test
    procedure index_r4                          ! compute interpolation index
    procedure inside_r4                         ! is it inside?
    procedure interpolate_r4                    ! interpolation
  end type
  integer:: id=0
CONTAINS

!> =============================================================================
!> Create a number of particles proportional at the left hand side wall
!> =============================================================================
SUBROUTINE init (self)
  class(pebble_patch_t):: self
  !.............................................................................
  character(len=8), save:: type='kdk'
  real:: density_yr_cgs = 1e-12, solid_density_cgs=3., size(3), pebble_size_cm=1.
  real:: t_max=1000.
  logical:: do_pebbles
  integer:: iostat
  namelist /pebble_params/ do_pebbles, type, density_yr_cgs, pebble_size_cm, &
    solid_density_cgs, size, t_max
  !-----------------------------------------------------------------------------
  ! Make sure to initialize, but only once
  !-----------------------------------------------------------------------------
  if (.not.associated(self%pebbles)) allocate (self%pebbles)
  if (self%pebbles%initialized) return
  self%pebbles%initialized = .true.
  call trace_begin ('pebbles_t%init', 3)
  call random%init
  call self%patch_t%init
  do_pebbles = .true.
  associate (mx=>self%mesh(1)%gn,my=>self%mesh(2)%gn,mz=>self%mesh(3)%gn,mt=>self%nt)
    allocate (self%density(mx,my,mz,mt))
    allocate (self%pressure(mx,my,mz,mt))
    allocate (self%vgas(mx,my,mz,3,mt))
  end associate
  print*,'shape',shape(self%density)
  print*,'shape',shape(self%pressure)
  print*,'shape',shape(self%vgas)
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
  self%pebbles%t_max = t_max
  self%pebbles%density = solid_density_cgs/scaling%d
  self%pebbles%mass = self%pebbles%density*self%pebbles%size**3
  self%pebbles%drag_coeff = scaling%d/(scaling%l*density_yr_cgs*cgs%yr/scaling%t)
  if (io%master) &
    print*,'pebbles: size, density, mass =', &
    self%pebbles%size, self%pebbles%density, self%pebbles%mass
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

!> =============================================================================
!> Create a number of particles proportional at the left hand side wall,
!> with initially a velocity interpolated from the argument velocity field
!> =============================================================================
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

!> =============================================================================
!> Kick step of kick-drift-kick (KDK), with gravity from star and planet
!> =============================================================================
SUBROUTINE kdk (self, p)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: p
  !.............................................................................
  real:: ar, rp(3)
  !-----------------------------------------------------------------------------
  rp = - p%p
  ar = 1./sqrt(sum(rp**2))
  p%ag = scaling%grav*planet%mass*rp*ar**3
  rp = self%pebbles%position - p%p
  ar = 1./sqrt(sum(rp**2))
  !p%ag = p%ag + scaling%grav*star%mass*rp*ar**3
  p%dt = self%pebbles%courant*sqrt(sum(p%v**2)/sum(p%ag**2))
  p%v = p%v + p%dt*p%ag

  p%p = p%p + 2.*p%dt*p%v
  p%t = p%t + 2.*p%dt

  rp = planet%position - p%p
  ar = 1./sqrt(sum(rp**2))
  p%ag = scaling%grav*planet%mass*rp*ar**3
  rp = star%position - p%p
  ar = 1./sqrt(sum(rp**2))
  !p%ag = p%ag + scaling%grav*star%mass*rp*ar**3
  p%v = p%v + p%dt*p%ag
END SUBROUTINE kdk

!> =============================================================================
!> Look up gas velocity components, given positions, and compute the drag rate
!> (inverse stopping time)
!> =============================================================================
SUBROUTINE drag (self, pebble)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: pebble
  !.............................................................................
  real:: w(3), d, p, temp, gamma, v(3), speed2, sound2
  integer:: ii(3), j, k
  integer, save:: nprint=100
  !-----------------------------------------------------------------------------
  call trace_begin ('pebbles_t%drag', 5)
  ii        = self%index_r4 (pebble%p, w, self%mesh%gn)
  d         = self%interpolate_r4 (self%density , ii, w)
  p         = self%interpolate_r4 (self%pressure, ii, w)
  if (any(w<0.0 .or. w>1.0)) &
    print*,'WARNING 1', ii, w
  !-----------------------------------------------------------------------------
  ! By using a 4-D pointer we can treat the velocity component as a 3D+time
  ! array, which can be passed to interpolate_r4.
  !-----------------------------------------------------------------------------
  block
    real, dimension(:,:,:,:), pointer:: velocity
    do j=1,3
      do k=1,3
        w(k) = (pebble%p(k) - self%mesh(k)%h(k+2))/self%mesh(k)%d + self%mesh(k)%o
        ii(k) = w(k)
        ii(k) = max(1,min(ii(k),self%mesh(k)%gn))
        w(k) = w(k) - ii(k)
      end do
      if (any(w<0.0 .or. w>1.0)) &
        print*,'WARNING 2', ii, w, j
      velocity => self%vgas(:,:,:,j,:)
      pebble%vg = self%interpolate_r4 (velocity, ii, w)
    end do
  end block
  speed2    = sum((pebble%v-pebble%vg)**2)
  temp      = p/d
  gamma     = merge(1.6666667,1.4,temp>2000.)
  sound2    = gamma*p/d
  pebble%c  = self%pebbles%drag_coeff*d/pebble%s*sqrt((speed2+sound2)/sound2)
  pebble%ad = (pebble%vg-pebble%v)*pebble%c
  if (nprint>0) then
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
  end if
  call trace_end
END SUBROUTINE drag

!> =============================================================================
!> Move a pebbles forward one step, with drag
!> =============================================================================
SUBROUTINE kdk_drag (self, p)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: p
  !.............................................................................
  real:: rp(3), ar
  !-----------------------------------------------------------------------------
  call trace_begin ('pebbles_t%kdk_drag', 5)
  !---------------------------------------------------------------------------
  ! Get acceleration from planet and star, ag(t)
  !---------------------------------------------------------------------------
  rp = - p%p
  ar = 1./sqrt(sum(rp**2))
  p%ag = scaling%grav*planet%mass*rp*ar**3
  rp = self%pebbles%position - p%p
  ar = 1./sqrt(sum(rp**2))
  p%ag = p%ag + scaling%grav*star%mass*rp*ar**3
  p%dt = self%courant*sqrt(sum(p%v**2)/sum(p%ag**2))

  !---------------------------------------------------------------------------
  ! Get the drag acceleration, at t+dt, then back up position
  !---------------------------------------------------------------------------
  p%p = p%p + p%dt*p%v
  p%t = p%t + p%dt
  call self%drag (p)
  p%p = p%p - p%dt*p%v

  !---------------------------------------------------------------------------
  ! First kick, with drag, v(t) -> v(t+dt)
  !---------------------------------------------------------------------------
  p%v = p%v + p%dt*(p%ag+p%c*p%dt*p%ad)/(1.+p%c*p%dt)

  !---------------------------------------------------------------------------
  ! Drift, r -> r+2*dt*v(t+dt*a), with v -> vg(t+dt) if drag dominates
  !---------------------------------------------------------------------------
  p%p = p%p + 2.*p%dt*p%v
  p%t = p%t + p%dt

  !---------------------------------------------------------------------------
  ! Get acceleration from planet and star, at p(t+2*dt)
  !---------------------------------------------------------------------------
  rp = - p%p
  ar = 1./sqrt(sum(rp**2))
  p%ag = scaling%grav*planet%mass*rp*ar**3
  rp = self%pebbles%position - p%p
  ar = 1./sqrt(sum(rp**2))
  p%ag = p%ag + scaling%grav*star%mass*rp*ar**3

  !---------------------------------------------------------------------------
  ! Get the drag acceleration, at t+2*dt
  !---------------------------------------------------------------------------
  call self%drag (p)

  !---------------------------------------------------------------------------
  ! Second kick, with drag, moving v to v(t+2*dt)
  !---------------------------------------------------------------------------
  p%v = p%v + p%dt*(p%ag+p%c*p%dt*p%ad)/(1.+p%c*p%dt)
  !---------------------------------------------------------------------------
  ! Next particle
  !---------------------------------------------------------------------------
  call trace_end
END SUBROUTINE kdk_drag

!> =============================================================================
!> Append a pebble to the list
!> =============================================================================
SUBROUTINE append (self, pebble)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble
  !.............................................................................
  pebble%prev => self%tail
  if (.not.associated(self%head)) self%head => pebble
  if (associated(self%tail)) self%tail%next => pebble
  self%tail => pebble
  self%n = self%n + 1
END SUBROUTINE append

!> =============================================================================
!> Remove a pebble from the list
!> =============================================================================
SUBROUTINE remove (self, pebble)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble
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

!> =============================================================================
SUBROUTINE export (self, pebble, nbor_list)
  class(pebble_patch_t):: self
  class(pebble_t), pointer:: pebble
  class(pebbles_t), pointer:: nbor_list
  !.............................................................................
  call self%pebbles%remove (pebble)
  call nbor_list%append (pebble)
END SUBROUTINE export

!> =============================================================================
!> Move the pebbles forward one sweep, and check if still inside our patch.
!> If not, search our neighbor list, and export to the neighbor pebble list.
!> Also check if it has hit the planet, and count that as an accretion
!> =============================================================================
SUBROUTINE update_pebbles (self, nbors)
  class(pebble_patch_t):: self
  class(link_t), pointer:: nbors
  !.............................................................................
  class(link_t), pointer:: nbor
  class(pebble_t), pointer:: pebble
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace_begin ('pebble_patch_t%update_pebbles',5)
  pebble => self%pebbles%head
  do while (associated(pebble))
    call self%move (pebble)
    if (.not. self%inside_r4 (pebble%p)) then
      if (associated(nbors)) then
        nbor => nbors
        do while (associated(nbor))
          task => nbor%task
          print*,'NBOR!'
          stop
          select type (task)
          class is (pebble_patch_t)
            if (task%inside_r4(pebble%p)) then
              call self%export (pebble, task%pebbles)
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
      end if
    else if (planet%inside(pebble%p)) then
      print*,'PLANET!', pebble%id
      call self%pebbles%remove (pebble)
      self%pebbles%n_accrete = self%pebbles%n_accrete + 1
    else if (pebble%t > self%pebbles%t_max) then
      call self%pebbles%remove (pebble)
    end if
    pebble => pebble%next
  end do
  self%pebbles%n_update = self%pebbles%n_update + self%pebbles%n
  call trace_end
END SUBROUTINE update_pebbles

!> =============================================================================
!> =============================================================================
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

!> =============================================================================
!> =============================================================================
FUNCTION inside_r4 (self, p) RESULT (out)
  class(pebble_patch_t):: self
  real, intent(in):: p(3)
  logical:: out
  !-----------------------------------------------------------------------------
  out = all(abs(p) < 0.5*self%size)
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

!> =============================================================================
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
!> =============================================================================
SUBROUTINE test (self, link)
  class(pebble_patch_t):: self
  !.............................................................................
  class(pebbles_t), pointer:: pebbles
  class(pebble_t), pointer:: pebble
  class(link_t), pointer:: link
  real:: position(3), size(3)
  integer:: k, it
  !-----------------------------------------------------------------------------
  call trace_begin ('pebbles_t%test')
  !-----------------------------------------------------------------------------
  ! Specify a star and a planet, and a 500 R_Earth box
  !-----------------------------------------------------------------------------
  call scaling%init
  star%mass   = cgs%m_sun/scaling%m
  planet%mass = cgs%m_earth/scaling%m
  star%position = 0.0
  planet%position = [0.0,1.0,0.0]*cgs%au/scaling%l
  planet%radius = cgs%r_earth/scaling%l
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
  !-----------------------------------------------------------------------------
  ! Gas properties
  !-----------------------------------------------------------------------------
  self%density = 1e-10/scaling%d
  self%pressure = self%density*10.
  print*,'TEMPERATURE',self%pressure(3,3,3,1)/self%density(3,3,3,1)*scaling%temp
  print*,shape(self%vgas)
  do it=1,self%nt
    self%vgas(:,:,:,1,it) = 0.5*cgs%kms/scaling%u*it
    self%vgas(:,:,:,2,it) = 0.0
    self%vgas(:,:,:,3,it) = 0.0
  end do
  block
    real, dimension(:,:,:,:), pointer:: velocity
    print*,associated(self%vgas)
    do it=1,self%nt
      velocity => self%vgas(:,:,:,1,:)
      print*,velocity(3,3,3,it)
      print*,associated(self%vgas)
      print*, 'vg_min,max', minval(self%vgas(:,:,:,1,it)), maxval(self%vgas(:,:,:,1,it))
      print*, 've_min,max', minval(velocity(:,:,:,it)), maxval(velocity(:,:,:,it))
    end do
  end block
  if (io%master) then
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
  call self%add (3000)
  self%move => kdk
  call tic
  it = 0
  do while (self%pebbles%n > 0)
    call self%update_pebbles (link%nbor)
    if (mod(it,10000)==0) then
      print *,self%pebbles%n, it
    end if
    it = it+1
  end do
  call toc ('kdk', self%pebbles%n_update)
  print*,'KDK accreted:', self%pebbles%n_accrete
  !-----------------------------------------------------------------------------
  ! Test the KDK + drag, with a 0.5 km/s head wind, blowing left to right
  !-----------------------------------------------------------------------------
  self%pebbles%n_update = 0
  self%pebbles%n_accrete = 0
  self%pebbles%t_max = 1e4
  call self%add (3000)
  self%move => kdk_drag
  do while (self%pebbles%n > 0)
    call self%update_pebbles (link%nbor)
  end do
  print*,'KDK_drag accreted', self%pebbles%n_accrete
  call trace_end
END SUBROUTINE test

END MODULE pebbles_mod
