!===============================================================================
!> Particles mover
!===============================================================================
MODULE particles_mod
  USE io_mod
  USE trace_mod
  USE omp_lock_mod
  USE timer_mod
  USE patch_mod
  USE link_mod
  USE kinds_mod
  implicit none
  private
  !.............................................................................
  ! EM-scaling; See Leung, Eur. J. Phys. 25 (2004)
  !.............................................................................
  type, public:: elm_t
    integer:: l(3), u(3)
    real:: kf=1., kB=1., kE=1.
  end type
  !.............................................................................
  ! nbor particles
  !.............................................................................
  type, public:: pnbor_t
    type(pnbor_t), pointer:: next => null()
    type(particles_t), pointer:: imports
  end type
  !.............................................................................
  ! Particles
  !.............................................................................
  type, public:: particles_t
    class(patch_t),    pointer :: patch   => null()
    type(particles_t), pointer :: exports => null()
    type(particles_t), pointer :: imports => null()
    type(pnbor_t),     pointer :: pnbors  => null()
    type(lock_t)               :: lock
    logical                    :: first_time = .true.
    integer                    :: id = 0
    integer                    :: m = 0
    integer                    :: n = 0
    integer                    :: l = 0
    integer                    :: offset(3) = 0
    integer                    :: seek = 0
    integer                    :: nprint = 0
    integer                    :: per_cell = 10
    integer                    :: lb(3), ub(3)
    real                       :: mass = 1.0
    real                       :: charge = 1.0
    real                       :: hB=1, hE=1, oc=1, oc2=1
    integer(2), pointer, contiguous :: q(:,:)         ! integer position
    real,       pointer, contiguous :: r(:,:)         ! float position
    real,       pointer, contiguous :: p(:,:)         ! momentum
    real,       pointer, contiguous :: e(:)           ! energy
    real,       pointer, contiguous :: w(:)           ! weight
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
    procedure:: update
    procedure:: post_update
    procedure, private:: init_light
    procedure, private:: export
    procedure, private:: import
    procedure, private:: move
  end type
  real    :: mass = 1.0
  real    :: charge = 1.0
  integer :: n_species = 2
  integer :: verbose = 1
  logical :: test = .false.
  type(elm_t), public:: elm
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE pre_init (self, patch, per_cell)
  class(particles_t):: self
  class(patch_t), pointer:: patch
  integer:: per_cell, iostat
  integer, save:: nprint=0
  namelist /particles_params/ mass, charge, test, nprint, verbose
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (self%first_time) then
    self%first_time = .false.
    read (io%input, particles_params, iostat=iostat)
    call io%namelist_warning ('particles_params', iostat)
    write (io%namelist, particles_params)
  end if
  call io%lock%unset
  !-----------------------------------------------------------------------------
  self%patch    => patch
  self%mass     = mass
  self%charge   = charge
  self%per_cell = per_cell
  self%lb       = patch%li
  self%ub       = patch%ui
  self%n        = per_cell*product(self%ub-self%lb+1)
  self%nprint   = nprint
END SUBROUTINE pre_init

!===============================================================================
!> Allocate particle arrays, with default values for tests
!===============================================================================
SUBROUTINE init (self, patch, n)
  class(particles_t):: self
  class(patch_t), pointer:: patch
  integer:: n
  !-----------------------------------------------------------------------------
  allocate (self%r(self%n,3))
  allocate (self%q(self%n,3))
  allocate (self%p(self%n,3))
  allocate (self%w(self%n))
  allocate (self%e(self%n))
  self%r = 0.0
  self%p = 0.0
  self%q(:,1) = self%lb(1)
  self%q(:,2) = self%lb(2)
  self%q(:,3) = self%lb(3)
  self%w = 1.0
END SUBROUTINE init

!===============================================================================
!> Import particles placed in our immports buffer by other tasks
!===============================================================================
SUBROUTINE init_light (self, patch, n)
  class(particles_t):: self
  class(patch_t), pointer:: patch
  integer:: n
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%init_light', itimer=itimer)
  call trace%end (itimer)
END SUBROUTINE init_light

!===============================================================================
!> Import particles and possible perform other preparations
!===============================================================================
SUBROUTINE pre_update (self)
  class(particles_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%pre_update')
  call self%import ()
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> Vay pusher from PPcode
!===============================================================================
SUBROUTINE update (self, B4, E4, J)
  class(particles_t):: self
  class(patch_t), pointer:: patch
  real(var_k), contiguous, dimension(:,:,:,:) :: B4, E4, J
  !.............................................................................
  real, dimension(3):: r, dr, p, p0, B, E, p1, t, ods, dq
  real:: c, o, oc, oc2, og, dt, odt, t2, s, ut, pp, hB, hE, p2, sqrt2
  integer:: i, q(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%update', itimer=itimer)
  patch => self%patch
  dt = merge (0d0, patch%dtime, test)
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! Update unrolled, picking up B and E via indirect addressing
  !-----------------------------------------------------------------------------
  ods = 1./patch%ds
  odt = 1./max(dt,1e-30)
  hE  = 0.5*dt*self%charge/self%mass              ! perfectly time centered
  hB  = hE*elm%kf                                 ! B prefactor
  oc  = 1. ! ...                                  ! 1./c
  oc2 = 1. !...                                   ! 1./c**2
  sqrt2 = sqrt(2.)
  if (verbose > 0 .and. self%nprint > 0) then
    call timer%tic
  end if
  do i=1,self%n
    !---------------------------------------------------------------------------
    ! Momentum update
    !---------------------------------------------------------------------------
    q = self%q(i,:)
    B = B4(q(1),q(2),q(3),:)
    E = E4(q(1),q(2),q(3),:)
    p = self%p(i,:)
    p0 = p
    o = 1./sqrt(1.+sum(p**2)*oc2) 
    t = hB*B
    t2 = sum(t**2)
    p1(1) = p(1) + 2.*hE*E(1) + o*(p(2)*t(3)-p(3)*t(2))   
    p1(2) = p(2) + 2.*hE*E(2) + o*(p(3)*t(1)-p(1)*t(3))   
    p1(3) = p(3) + 2.*hE*E(3) + o*(p(1)*t(2)-p(2)*t(1))   
    s = 1. + sum(p**2)*oc2 - t2
    pp = sum(p1*t)*oc
    o = sqrt2 * 1./sqrt(s + sqrt(s*s + 4.*(t2 + pp*pp))) 
    t= t*o
    s = 1./(1. + t2*o*o)
    ut = sum(p1*t)
    p(1) = s * (p1(1) + ut*t(1) + p1(2)*t(3) - p1(3)*t(2))
    p(2) = s * (p1(2) + ut*t(2) + p1(3)*t(1) - p1(1)*t(3))
    p(3) = s * (p1(3) + ut*t(3) + p1(1)*t(2) - p1(2)*t(1))
    !---------------------------------------------------------------------------
    ! Energy update
    !---------------------------------------------------------------------------
    p2= sum(p**2)                                             ! momentum squared
    og = 1./sqrt(1.+p2*oc2)                                   ! 1./gamma
    self%e(i) = self%mass*p2*og/(1.+og)                       ! kinetic energy
    !---------------------------------------------------------------------------
    ! Position update
    !---------------------------------------------------------------------------
    og = og*dt
    dr = p*og*ods
    r = self%r(i,:) + dr
    dq = floor(r)
    self%q(i,:) = q + dq
    self%r(i,:) = r - dq
    self%p(i,:) = p
    J(:,q(1),q(2),q(3)) = &
    J(:,q(1),q(2),q(3)) + self%charge*self%w(i)*dr*odt
  end do
  if (verbose > 0 .and. self%nprint > 0) then
    self%nprint = self%nprint-1
    call timer%toc ('Vay push', self%n)
  end if
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Export and possible renormalize particle weights, etc
!===============================================================================
SUBROUTINE post_update (self)
  class(particles_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%post_update')
  call self%export ()
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Export particles that left our ROA into imports buffers of nbor patches
!===============================================================================
SUBROUTINE export (self)
  class(particles_t):: self
  class(patch_t), pointer:: patch, nbpatch
  class(pnbor_t), pointer:: pnbor
  class(particles_t), pointer:: to, fr, exports
  class(nbor_t), pointer:: nbor
  integer:: i, j, offset(3), nbox(3), l(3), u(3), q(3), n 
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%exports%n < 1) &
    return
  call trace%begin ('particles_t%export', itimer=itimer)
  patch => self%patch
  exports => self%exports
if (verbose>0) PRINT *,exports%n,' particles in exports',exports%id
  call exports%init_light (patch, exports%n)
  exports%id = self%id
  self%imports%id = self%id
  !-----------------------------------------------------------------------------
  ! Move particles that fall outside our region of authority to exports
  !-----------------------------------------------------------------------------
  l = self%lb
  u = self%ub
  n = 0
  do i = 1,self%m
    if (self%w(i) < 0.0) then
      n = n+1
if (verbose>0) PRINT *,i,self%w(i),self%q(i,:)
    end if
  end do
if (verbose>0) PRINT *,n,' particles with negative w',exports%id
  n = 0
  j = 0
  do i = 1,self%m
    if (self%w(i) < 0.0) then
      self%w(i) = -self%w(i)
      n = n+1
      j = i
      call self%move (exports, i)
    end if
  end do
if (verbose>0) PRINT *,exports%n,' particles in exports',exports%id, n, j
  !-----------------------------------------------------------------------------
  ! Loop over nbor patches and copy exports particles that fall inside an nbor
  ! patch to its imports buffer
  !-----------------------------------------------------------------------------
  nbox = nint(patch%box/patch%ds)
  pnbor => self%pnbors
  do while(associated(pnbor))    
    to => pnbor%imports
    offset = exports%offset - to%offset
                                                     call to%lock%set ('export')
    do i=1,exports%n
      q = exports%q(i,:)
      q = modulo (q + offset, nbox)
      if (all(q >= l .and. q <= u)) then
        call exports%move (to, i)
      end if
    end do
if (verbose>0) PRINT *,to%n,' particles in imports to', to%id
                                                   call to%lock%unset ('export')
    pnbor => pnbor%next
  end do
if (verbose>0) PRINT *,exports%n,' particles left in exports before BCs'
  call trace%end (itimer)
END SUBROUTINE export

!===============================================================================
!> Import particles placed in our immports buffer by other tasks
!===============================================================================
SUBROUTINE import (self)
  class(particles_t):: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%import', itimer=itimer)
  call trace%end (itimer)
END SUBROUTINE import

!==============================================================================
!> Copy particle i from self to recv and increment n, if the slot is empty
!===============================================================================
SUBROUTINE move (self, recv, j)
  class(particles_t):: self, recv
  integer:: j, i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%n < 1) return
  call trace%begin ('particles_t%move', itimer=itimer)
  i = recv%seek
  do while (recv%w(i) > 0.0)
    i = i + 1    
  end do
  if (i <= recv%m) then
    recv%r(i,:) = self%r(j,:)
    recv%q(i,:) = self%q(j,:)
    recv%p(i,:) = self%p(j,:)
    recv%w(i)   = self%w(j)
    self%w(j)   = 0
    self%q(j,:) = self%lb
    recv%seek = i
    recv%l = max(i,recv%l)
    recv%n = recv%n + 1
    self%n = self%n - 1
  else
    write (stderr,*) 'particles_t%move: out of space'
  end if
  call trace%end (itimer)
END SUBROUTINE move

END MODULE particles_mod
