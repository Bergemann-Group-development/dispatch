!===============================================================================
!> MUSCL method using Riemann solver to compute fluxes
!>
!> The Riemann solver takes pointers to the left and right interface values
!> as input, derived from the primitive variables and their slopes.
!>
!> As primitive variables we choose log density, internal energy per unit mass,
!> and velocities. This ensures that the densities are positive, and that the
!> slopes are derived from well-conditioned per-unit-mass variables.
!===============================================================================
MODULE solver_mod
  USE io_mod
  USE riemann_mod
  USE omp_lib
  implicit none
  private
  integer, parameter:: nt=5
  type, public:: solver_t
    integer:: id
    integer:: nw=1, nv=5, nt=nt, ng=2, n=16
    integer:: it=1, new=2
    integer, dimension(3):: gn, lb, lo, li, ui, uo, ub
    real(8):: time=0d0, dtime=1d0
    real(8):: t(nt), dt(nt), ds(3)=1d0
    real:: slope_type
    real:: courant
    real:: csound
    real(8):: g1, g2, gamma
    character(len=16):: kind = 'ramses_hd_patch'
    !type(riemann_t):: riemann
    real, dimension(:,:,:,:,:,:),  allocatable:: mem
    real, dimension(:,:,:,:)  , allocatable :: prim
    real, dimension(:,:,:,:)  , allocatable        :: left
    real, dimension(:,:,:,:)  , allocatable        :: rght
    real, dimension(:,:,:,:,:), allocatable :: grad
    real, dimension(:,:,:,:,:), allocatable        :: flux
  contains
    procedure:: init
    procedure:: update
  end type
  real(8):: start
  !-----------------------------------------------------------------------------
  ! Input parameters:
  !-----------------------------------------------------------------------------
  integer:: n=16                        ! internal cells
  integer:: ng=2                        ! guard cells
  integer:: slope_type=2                ! slope limiter
  integer:: verbose=0                   ! info level
  logical:: detailed_timer=.true.       ! optional timer
  real(8):: gamma=1.4                   ! gas gamma
  real:: csound=1.0                     ! sound speed
  real:: courant=0.2                    ! courant number
  character(len=4):: method='hll'       ! riemann method

  !$omp declare target(csound, slope_type, method)
  !$omp declare target(ctoprim, sources, divflux, riemn3d, predict, slopes)
CONTAINS

!===============================================================================
!> Initialize parameters
!===============================================================================
SUBROUTINE init (self)
  class(solver_t):: self
  !.............................................................................
  namelist /solver_params/ verbose, courant, gamma, csound, slope_type, &
    method, detailed_timer
  integer:: iostat, i
  integer:: m
  integer, save:: id=0
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  !$omp atomic update
  id = id + 1
  self%id = id
  !-----------------------------------------------------------------------------
  ! Here we replace the call patch_t%init directly, to initialize %mem
  !-----------------------------------------------------------------------------
  self%gn = self%n + 2*ng    ! all cells
  self%lb = 1                ! lower buondary
  self%li = 1 + ng           ! lower internal index
  self%lo = self%li - 1      ! lower outside index
  self%ui = self%lo + n      ! upper internal index
  self%uo = self%ui + 1      ! upper outside index
  self%ub = self%ui + ng     ! uppper boundary
  m = self%gn(1)
  allocate (self%mem(m,m,m,self%nv,self%nt,self%nw))
  self%mem(:,:,:,1,:,:) = 1.0
  self%mem(:,:,:,2,:,:) = 0.0
  self%mem(:,:,:,3,:,:) = 0.0
  self%mem(:,:,:,4,:,:) = 0.0
  self%mem(:,:,:,5,:,:) = 1.0
  call io%lock%set
  if (first_time) then
    print *, 'solver_t%init:'
    print *, 'id', self%id
    print *, 'gn', self%gn
    print *, 'lb', self%lb
    print *, 'lo', self%lo
    print *, 'li', self%li
    print *, 'ub', self%ub
    print *, 'uo', self%uo
    print *, 'ui', self%ui
    print *, 'mem', shape(self%mem)
    print *, 'muscl2'
    start = omp_get_wtime()
    first_time = .false.
  end if
  call io%lock%unset
  !call self%patch_t%init
  !-----------------------------------------------------------------------------
  ! Ramses-like specific things
  !-----------------------------------------------------------------------------
  !self%riemann%method = method
  self%courant = courant
  self%csound = csound
  !-----------------------------------------------------------------------------
  ! Make sure to save and use consistent gamma factors throughout
  !-----------------------------------------------------------------------------
  self%gamma = gamma
  if (gamma == 1d0) then
    self%g1 = 1d0
  else
    self%g1 = gamma-1d0
  end if
  self%g2 = self%g1*gamma
  !self%riemann%g1 = self%g1
  !self%riemann%g2 = self%g2
  call alloc   (self)
  !-----------------------------------------------------------------------------
  ! Generic settings
  !-----------------------------------------------------------------------------
  !self%staggered = .false.
  !do i=1,3
  !  self%mesh(i)%h = 0d0
  !end do
  !self%unsigned(self%idx%d) = unsigned
  !self%unsigned(self%idx%e) = unsigned
  !self%pervolume(self%idx%px) = unsigned
  !self%pervolume(self%idx%py) = unsigned
  !self%pervolume(self%idx%pz) = unsigned
  !-----------------------------------------------------------------------------
  ! Calling gpatch_t%init, which may call ICs, must wait until all properties
  ! of patches (mesh, idx, patch_t%h staggering) are in place.
  !-----------------------------------------------------------------------------
  !call self%gpatch_t%init
END SUBROUTINE init

!===============================================================================
!> Allocate work arrays
!===============================================================================
SUBROUTINE update (self)
  class(solver_t), target:: self
  integer, save:: itimer=0
  integer, dimension(3) :: lb, ub, lo, uo, li, ui
  real, dimension(:,:,:,:,:), pointer :: grad
  real, dimension(:,:,:,:), pointer :: prim
  real, dimension(:,:,:,:,:,:),  pointer:: mem
  real, dimension(:,:,:,:,:), pointer :: flux
  real, dimension(:,:,:,:)  , pointer :: left
  real, dimension(:,:,:,:)  , pointer :: rght
  integer :: new, it, nv
  real(8) :: gamma, ds(3), dtime, g1,g2
  !type(riemann_t) :: riemann
  !.............................................................................
  !call self%output
  !call alloc   (self)
  grad => self%grad
  prim => self%prim
  mem => self%mem
  flux => self%flux
  left => self%left
  rght => self%rght
  nv = self%nv
  lb = self%lb
  ub = self%ub
  new = self%new
  it = self%it
  gamma = self%gamma
  lo = self%lo
  uo = self%uo
  dtime = self%dtime
  ds = self%ds
  g1 = self%g1
  g2 = self%g2
  li = self%li
  ui = self%ui
  !riemann = self%riemann
  !$omp target data map(to:grad) map(to:prim) map(tofrom:mem) map(to: flux) map(to:left) map(to:rght) map(to:self)
  !$omp target teams
  call ctoprim (self, mem, new, it, prim, lb, ub, gamma)
  call slopes  (self, nv, prim, grad, lb, ub)
  call predict (self, grad, prim, mem, new, dtime, ds, g1, lo, uo)
  call riemn3d (self, prim, grad, flux, left, rght, nv, lo, uo, g1, g2)
  call divflux (self, mem, flux, nv, new, ds, dtime, li, ui)
  call sources (self, mem, prim, new, gamma)
  !$omp end target teams
  !$omp end target data
  !call dealloc (self)
  self%time = self%time + self%dtime
  print '(f12.6,2x,a,i4,f8.2)', &
    omp_get_wtime()-start, 'solver_t%update: id, time =', self%id, self%time
  !call self%counter_update
END SUBROUTINE update

!===============================================================================
!> Allocate work arrays
!===============================================================================
SUBROUTINE alloc (self)
  class(solver_t):: self
  integer, save:: itimer=0
  !.............................................................................
  call tracebegin ('solver_t%alloc', itimer=itimer)
  if (.not.allocated(self%prim)) then

    allocate (self%prim (self%gn(1), self%gn(2), self%gn(3), self%nv   ))
    allocate (self%left (self%gn(1), self%gn(2) ,self%gn(3), self%nv   ))
    allocate (self%rght (self%gn(1), self%gn(2) ,self%gn(3), self%nv   ))
    allocate (self%grad (self%gn(1), self%gn(2) ,self%gn(3), self%nv, 3))
    allocate (self%flux (self%gn(1), self%gn(2) ,self%gn(3), self%nv, 3))
    !---------------------------------------------------------------------------
    ! Set to zero once, to allow undefined guard zone values to be used
    !---------------------------------------------------------------------------

    self%prim = 0.0
    self%left = 0.0
    self%rght = 0.0
    self%grad = 0.0
    self%flux = 0.0
  end if
  call traceend (itimer)
END SUBROUTINE alloc

!===============================================================================
!> Deallocate work arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(solver_t):: self
  integer, save:: itimer=0
  !.............................................................................
  deallocate (self%prim)
  deallocate (self%left)
  deallocate (self%rght)
  deallocate (self%grad)
  deallocate (self%flux)
END SUBROUTINE dealloc

!===============================================================================
!> Conservative to primitive transformation, using log(d), U, and Eint
!> Valid in full interval: [lb,ub] = [li-2,ui+2]
!===============================================================================
SUBROUTINE ctoprim (self, mem, new, it, prim, l, u, gamma)
  class(solver_t), target:: self
  !.............................................................................
  real, dimension(:,:,:,:)  , pointer :: prim
  real, dimension(:,:,:,:,:,:), pointer :: mem
  integer:: i1, i2, i3, iv, l(3), u(3)
  !type(lock_t), pointer:: lock
  real:: u2, u2_max, c2, g2
  integer, save:: itimer=0
  integer :: new
  integer :: it
  real(8) :: gamma
  !new = self%new
  !prim => self%prim
  !mem => self%mem
  !-----------------------------------------------------------------------------
  !call tracebegin ('solver_t%ctoprim', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Copy everything to new time slot, lock protected
  !-----------------------------------------------------------------------------
  !lock => self%mem_lock(self%it)
  !call lock%set ('ctoprim')
  !$omp parallel workshare
  mem(:,:,:,:,new,1) = mem(:,:,:,:,it,1)
  !$omp end parallel workshare
  !call lock%unset ('ctoprim')
  !-----------------------------------------------------------------------------
  !l = self%lb
  !u = self%ub
  !$omp distribute parallel do simd collapse(2)
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    prim(:,i2,i3,1) = log(mem(:,i2,i3, 1,new,1))
    do iv=2,5
      prim(:,i2,i3,iv)  = mem(:,i2,i3,iv,new,1) &
                        / mem(:,i2,i3, 1,new,1)
    end do
  end do
  end do
  !$omp end distribute parallel do simd
  !-----------------------------------------------------------------------------
  ! Subtract kinetic energy per unit mass and compute self%u_max
  !-----------------------------------------------------------------------------
  if (gamma == 1d0) then
    g2 = 1.0
  else
    g2 = (gamma-1d0)*gamma
  end if
  u2_max = 0.0
  !$omp distribute parallel do simd collapse(3)
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    do i1=l(1),u(1)
      u2 = prim(i1,i2,i3,2)**2 + &
           prim(i1,i2,i3,3)**2 + &
           prim(i1,i2,i3,4)**2
      prim(i1,i2,i3,5) = prim(i1,i2,i3,5) - 0.5*u2
      c2 = prim(i1,i2,i3,5)*g2
      u2_max = max (u2_max, u2+c2)
      prim(i1,i2,i3,5) = log(prim(i1,i2,i3,5))
    end do
  end do
  end do
  !$omp end distribute parallel do simd
  !self%u_max = sqrt(u2_max)
!print *,self%id,self%u_max
  !call self%patch_t%courant_condition
  !call traceend (itimer)
END SUBROUTINE ctoprim

!===============================================================================
!> Slope limiters, produce slopes in the interval [lb+1,ub-1] = [lo,uo]
!>
!> Costs at 16**3 point, in ms/call (divide with 4 to get mus/cell):
!>            T430
!> central    0.08
!> minmod     0.24
!> moncen     0.37
!===============================================================================
SUBROUTINE slopes (self, nv, prim, grad, lb, ub)
  class(solver_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, dir, l(3), u(3), lb(3), ub(3)
  integer, save:: itimer=0
  integer :: nv
  real, dimension(:,:,:,:)  , pointer :: prim
  real, dimension(:,:,:,:,:)  , pointer :: grad
  !$omp declare target(central, minmod, moncen)
  !-----------------------------------------------------------------------------
  !call tracebegin ('solver_t%slopes', itimer=itimer)
  l = lb+1
  u = ub-1
  if (slope_type == 1) then
    call minmod (self)
  else if (slope_type >= 2) then
    call moncen (self)
  else
    call central (self)
  end if
  !call traceend (itimer)
contains
!===============================================================================
subroutine central (self)
  class(solver_t), target:: self
  integer :: nv
  real, dimension(:,:,:,:)  , pointer :: prim
  real, dimension(:,:,:,:,:)  , pointer :: grad

  prim => self%prim
  grad => self%grad
  nv = self%nv
  !$omp distribute parallel do collapse(4)
  do iv=1,nv
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    do i1=l(1),u(1)
      grad (i1,i2,i3,iv,1) = 0.5*(prim(i1+1,i2,i3,iv) - prim(i1-1,i2,i3,iv))
      grad (i1,i2,i3,iv,2) = 0.5*(prim(i1,i2+1,i3,iv) - prim(i1,i2-1,i3,iv))
      grad (i1,i2,i3,iv,3) = 0.5*(prim(i1,i2,i3+1,iv) - prim(i1,i2,i3-1,iv))
    end do
  end do
  end do
  end do
  !$omp end distribute parallel do
end subroutine central
!===============================================================================
subroutine minmod (self)
  class(solver_t), target:: self
  real:: dlft, drgt
  integer :: nv
  real, dimension(:,:,:,:)  , pointer :: prim
  real, dimension(:,:,:,:,:)  , pointer :: grad

  prim => self%prim
  grad => self%grad
  nv = self%nv

  !$omp distribute parallel do collapse(4)
  do iv=1,nv
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    do i1=l(1),u(1)
      ! slopes in first coordinate direction
      dlft = prim(i1  ,i2,i3,iv) - prim(i1-1,i2,i3,iv)
      drgt = prim(i1+1,i2,i3,iv) - prim(i1  ,i2,i3,iv)
      grad(i1,i2,i3,iv,1) = merge(0.0,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=0.0)
      ! slopes in second coordinate direction
      dlft = prim(i1,i2  ,i3,iv) - prim(i1,i2-1,i3,iv)
      drgt = prim(i1,i2+1,i3,iv) - prim(i1,i2  ,i3,iv)
      grad(i1,i2,i3,iv,2) = merge(0.0,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=0.0)
      ! slopes in third coordinate direction
      dlft = prim(i1,i2,i3  ,iv) - prim(i1,i2,i3-1,iv)
      drgt = prim(i1,i2,i3+1,iv) - prim(i1,i2,i3  ,iv)
      grad(i1,i2,i3,iv,3) = merge(0.0,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=0.0)
    end do
  end do
  end do
  end do
  !$omp end distribute parallel do
end subroutine minmod
!===============================================================================
subroutine moncen(self)
  class(solver_t), target:: self
  real:: dlft,  drgt, dcen, dsgn, slop, dlim, c1, c2
  !integer :: nv
  !real, dimension(:,:,:,:)  , pointer :: prim
  !real, dimension(:,:,:,:,:)  , pointer :: grad

  !prim => self%prim
  !grad => self%grad
  !nv = self%nv
  c2 = 0.5
  !$omp distribute parallel do collapse(3)
  do iv=1,nv
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !$omp simd private(dlft,drgt,dcen,dsgn,slop,dlim)
    do i1 = l(1),u(1)
      ! slopes in first coordinate direction
      dlft = prim(i1  ,i2,i3,iv) - prim(i1-1,i2,i3,iv)
      drgt = prim(i1+1,i2,i3,iv) - prim(i1  ,i2,i3,iv)
      dcen = c2*(dlft+drgt)
      dsgn = sign(1.0, dcen)
      slop = min(abs(dlft),abs(drgt))
      dlim = merge(0.0,slop,(dlft*drgt)<=0.0)
      grad(i1,i2,i3,iv,1) = dsgn*min(dlim,abs(dcen))
      ! slopes in second coordinate direction
      dlft = prim(i1,i2  ,i3,iv) - prim(i1,i2-1,i3,iv)
      drgt = prim(i1,i2+1,i3,iv) - prim(i1,i2  ,i3,iv)
      dcen = c2*(dlft+drgt)
      dsgn = sign(1.0,dcen)
      slop = min(abs(dlft),abs(drgt))
      dlim = merge(0.0,slop,(dlft*drgt)<=0.0)
      grad(i1,i2,i3,iv,2) = dsgn*min(dlim,abs(dcen))
      ! slopes in third coordinate direction
      dlft = prim(i1,i2,i3  ,iv) - prim(i1,i2,i3-1,iv)
      drgt = prim(i1,i2,i3+1,iv) - prim(i1,i2,i3  ,iv)
      dcen = c2*(dlft+drgt)
      dsgn = sign(1.0,dcen)
      slop = min(abs(dlft),abs(drgt))
      dlim = merge(0.0,slop,(dlft*drgt)<=0.0)
      grad(i1,i2,i3,iv,3) = dsgn*min(dlim,abs(dcen))
    end do
  end do
  end do
  end do
  !$omp end distribute parallel do
end subroutine moncen
END SUBROUTINE slopes

!===============================================================================
!> Predictor step for self%primitive variables; should produce self%left and right states
!> in intervals [li,uo] for self%left and [lo,uo] for right, so that after a shift
!> of right one index, the right state is valid in the interval [lo,ui]
!===============================================================================
SUBROUTINE predict (self, grad, prim, mem, new, dtime, ds, g1, l, u)
  class(solver_t), target :: self
  real(8) :: dtime, ds(3), g1
  integer :: new
  real, dimension(:,:,:,:,:), pointer :: grad
  real, dimension(:,:,:,:), pointer :: prim
  real, dimension(:,:,:,:,:,:),  pointer:: mem
  !real, pointer, dimension(:,:,:):: &
  !  ln_1, ln_2, ln_3, &
  !  ee_1, ee_2, ee_3, &
  !  u1_1, u1_2, u1_3, &
  !  u2_1, u2_2, u2_3, &
  !  u3_1, u3_2, u3_3, &
  !  ln, ee, u1, u2, u3, d
  real:: divu, dtd1, dtd2, dtd3, dt2
  real:: ff(size(prim,1),3)
  integer:: i1, i2, i3, l(3), u(3)
  integer, save:: itimer=0
  !.............................................................................
  !call tracebegin ('solver_t%predict', itimer=itimer)
  !ln_1 => grad(:,:,:,1,1)
  !ln_2 => grad(:,:,:,1,2)
  !ln_3 => grad(:,:,:,1,3)
  !u1_1 => grad(:,:,:,2,1)
  !u1_2 => grad(:,:,:,2,2)
  !u1_3 => grad(:,:,:,2,3)
  !u2_1 => grad(:,:,:,3,1)
  !u2_2 => grad(:,:,:,3,2)
  !u2_3 => grad(:,:,:,3,3)
  !u3_1 => grad(:,:,:,4,1)
  !u3_2 => grad(:,:,:,4,2)
  !u3_3 => grad(:,:,:,4,3)
  !ee_1 => grad(:,:,:,5,1)
  !ee_2 => grad(:,:,:,5,2)
  !ee_3 => grad(:,:,:,5,3)
  !ln   => prim(:,:,:,1)
  !u1   => prim(:,:,:,2)
  !u2   => prim(:,:,:,3)
  !u3   => prim(:,:,:,4)
  !ee   => prim(:,:,:,5)
  !d    => mem(:,:,:,1,new,1)
  !-----------------------------------------------------------------------------
  ! Evolution of self%primitive variables over half a time step, with
  ! \partial\ln(d) / \partial t = -u*self%grad(ln_d) - div(u)
  !     \partial u / \partial t = -u*self%grad(u) - self%grad(P)/d
  !                     self%grad(P) = d*self%grad(e) + e*self%grad(d)
  !                             = d*self%grad(e) + e*d*self%grad(\ln(d))
  !     \partial e / \partial t = -u*self%grad(e) - (P/d)*div(U)
  !-----------------------------------------------------------------------------
  dt2  = dtime*0.5
  dtd1 = dt2/ds(1)
  dtd2 = dt2/ds(2)
  dtd3 = dt2/ds(3)
  !g1 = g1
  !l = lo
  !u = uo
  !$omp distribute parallel do simd collapse(3)
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !ff(:,:) = 0.0
    !if (allocated(self%force_per_unit_mass)) then
    !  ff(:,1) = ff(:,1) + self%force_per_unit_mass(:,i2,i3,1)
    !  ff(:,2) = ff(:,2) + self%force_per_unit_mass(:,i2,i3,2)
    !  ff(:,3) = ff(:,3) + self%force_per_unit_mass(:,i2,i3,3)
    !end if
    !if (allocated(self%force_per_unit_volume)) then
    !  ff(:,1) = ff(:,1) + self%force_per_unit_volume(:,i2,i3,1)/d(:,i2,i3)
    !  ff(:,2) = ff(:,2) + self%force_per_unit_volume(:,i2,i3,2)/d(:,i2,i3)
    !  ff(:,3) = ff(:,3) + self%force_per_unit_volume(:,i2,i3,3)/d(:,i2,i3)
    !end if
    !!$omp simd
    do i1=l(1),u(1)
      divu = grad(i1,i2,i3,2,1)*dtd1 + grad(i1,i2,i3,3,2)*dtd2 + grad(i1,i2,i3,4,3)*dtd3
      !divu = u1_1(i1,i2,i3)*dtd1 + u2_2(i1,i2,i3)*dtd2 + u3_3(i1,i2,i3)*dtd3
      prim(i1,i2,i3,1) = prim(i1,i2,i3,1) - ( &
        prim(i1,i2,i3,2)*grad(i1,i2,i3,1,1)*dtd1 + &
        prim(i1,i2,i3,3)*grad(i1,i2,i3,1,2)*dtd2 + &
        prim(i1,i2,i3,4)*grad(i1,i2,i3,1,3)*dtd3) - divu
      prim(i1,i2,i3,2) = prim(i1,i2,i3,2) - ( & !u1(i1,i2,i3) - ( &
        prim(i1,i2,i3,2)*grad(i1,i2,i3,2,1)*dtd1 + & !u1(i1,i2,i3)*u1_1(i1,i2,i3)*dtd1 + &
        prim(i1,i2,i3,3)*grad(i1,i2,i3,2,2)*dtd2 + & !u2(i1,i2,i3)*u1_2(i1,i2,i3)*dtd2 + &
        prim(i1,i2,i3,4)*grad(i1,i2,i3,2,3)*dtd3) - & !u3(i1,i2,i3)*u1_3(i1,i2,i3)*dtd3) - &
        (grad(i1,i2,i3,5,1) + prim(i1,i2,i3,5)*grad(i1,i2,i3,1,1))*g1*dtd1! + ff(i1,1)*dt2
        !(ee_1(i1,i2,i3) + ee(i1,i2,i3)*ln_1(i1,i2,i3))*g1*dtd1! + ff(i1,1)*dt2
      prim(i1,i2,i3,3) = prim(i1,i2,i3,3) - ( & !u2(i1,i2,i3) - ( &
        prim(i1,i2,i3,2)*grad(i1,i2,i3,3,1)*dtd1 + &  !u1(i1,i2,i3)*u2_1(i1,i2,i3)*dtd1 + &
        prim(i1,i2,i3,3)*grad(i1,i2,i3,3,2)*dtd2 + &  !u2(i1,i2,i3)*u2_2(i1,i2,i3)*dtd2 + &
        prim(i1,i2,i3,4)*grad(i1,i2,i3,3,3)*dtd3) - & !u3(i1,i2,i3)*u2_3(i1,i2,i3)*dtd3) - &
        (grad(i1,i2,i3,5,2) + prim(i1,i2,i3,5)*grad(i1,i2,i3,1,2))*g1*dtd2! + ff(i1,2)*dt2
        !(ee_2(i1,i2,i3) + ee(i1,i2,i3)*ln_2(i1,i2,i3))*g1*dtd2! + ff(i1,2)*dt2
      prim(i1,i2,i3,4) = prim(i1,i2,i3,4) - ( & !u3(i1,i2,i3) - ( &
        prim(i1,i2,i3,2)*grad(i1,i2,i3,4,1)*dtd1 + &  !u1(i1,i2,i3)*u3_1(i1,i2,i3)*dtd1 + &
        prim(i1,i2,i3,3)*grad(i1,i2,i3,4,2)*dtd2 + &  !u2(i1,i2,i3)*u3_2(i1,i2,i3)*dtd2 + &
        prim(i1,i2,i3,4)*grad(i1,i2,i3,4,3)*dtd3) - & !u3(i1,i2,i3)*u3_3(i1,i2,i3)*dtd3) - &
        (grad(i1,i2,i3,5,3) + prim(i1,i2,i3,5)*grad(i1,i2,i3,1,3))*g1*dtd3! + ff(i1,3)*dt2
        !(ee_3(i1,i2,i3) + ee(i1,i2,i3)*ln_3(i1,i2,i3))*g1*dtd3! + ff(i1,3)*dt2
      prim(i1,i2,i3,5) = prim(i1,i2,i3,5) - ( & !ee(i1,i2,i3) - ( &
        prim(i1,i2,i3,2)*grad(i1,i2,i3,5,1)*dtd1 + &        !u1(i1,i2,i3)*ee_1(i1,i2,i3)*dtd1 + &
        prim(i1,i2,i3,3)*grad(i1,i2,i3,5,2)*dtd2 + &        !u2(i1,i2,i3)*ee_2(i1,i2,i3)*dtd2 + &
        prim(i1,i2,i3,4)*grad(i1,i2,i3,5,3)*dtd3) - divu*g1 !u3(i1,i2,i3)*ee_3(i1,i2,i3)*dtd3) - divu*g1
    end do
  end do
  end do
  !$omp end distribute parallel do simd
  !call traceend (itimer)
END SUBROUTINE predict

!===============================================================================
!> Select the self%left and right states and call the Riemann solver to get self%fluxes.
!> The self%left and right states are defined for indices [lo,uo], which is where
!> the slopes are valid.  The self%left states refer to positions [lo+0.5,uo+0.5]
!> and the right states refer to positions [lo-0.5,uo-0.5], while we need the
!> self%fluxes in the interval [li-0.5,ui+0.5] = [lo+0.5,uo-0.5], which matches
!> what the self%left and right states combined can supply.  In order to refer to
!> self%left and right states at the same physical location, we stay with the
!> indices of the right states, whose physical locations match those where we
!> want the self%fluxes, while referring to the self%left state one index to the self%left
!> inside the Riemann solver
!===============================================================================
SUBROUTINE riemn3d (self, prim, grad, flux, left, rght, nv, l, u, g1, g2)
  !!$omp declare target
  class(solver_t), target:: self
  !.............................................................................
  integer:: l(3), u(3), i1, i2, i3, iv, nv
  integer, save:: itimer=0, itimer1=0
  real, dimension(:,:,:,:), pointer :: prim
  real, dimension(:,:,:,:,:), pointer :: grad
  real, dimension(:,:,:,:,:), pointer :: flux
  real, dimension(:,:,:,:)  , pointer :: left
  real, dimension(:,:,:,:)  , pointer :: rght
  real, dimension(:,:,:,:)  , pointer :: solverflux
  real(8):: g1,g2
  !type(riemann_t) :: riemann
  !-----------------------------------------------------------------------------
  !call tracebegin ('solver_t%riemn3d', itimer=itimer)
  !grad => self%grad
  !prim => self%prim
  !flux => self%flux
  !left => self%left
  !rght => self%rght
  !nv = self%nv
  !l = self%lo
  !u = self%uo
  !riemann = self%riemann
  !-----------------------------------------------------------------------------
  ! X-direction
  !-----------------------------------------------------------------------------
  !$omp distribute parallel do simd collapse(4)
  do iv=1,nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    left(i1,i2,i3,iv) = prim(i1,i2,i3,iv) + grad(i1,i2,i3,iv,1)*0.5
    rght(i1,i2,i3,iv) = prim(i1,i2,i3,iv) - grad(i1,i2,i3,iv,1)*0.5
  end do; end do; end do; end do
  !$omp end distribute parallel do simd
  !                  call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  !solverflux => flux(:,:,:,:,1)
  call solve (1, l+1, u, left, rght, flux,g1,g2)
  !                  call traceend (itimer1)
  !-----------------------------------------------------------------------------
  ! Y-direction
  !-----------------------------------------------------------------------------
  !$omp distribute parallel do simd collapse(4)
  do iv=1,nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    left(i1,i2,i3,iv) = prim(i1,i2,i3,iv) + grad(i1,i2,i3,iv,2)*0.5
    rght(i1,i2,i3,iv) = prim(i1,i2,i3,iv) - grad(i1,i2,i3,iv,2)*0.5
  end do; end do; end do; end do
  !$omp end distribute parallel do simd
  !                  call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  !solverflux => flux(:,:,:,:,2)
  call solve (2, l+1, u, left, rght, flux,g1,g2)
  !                  call traceend (itimer1)
  !-----------------------------------------------------------------------------
  ! Z-direction
  !-----------------------------------------------------------------------------
  !$omp distribute parallel do simd collapse(4)
  do iv=1,nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    left(i1,i2,i3,iv) = prim(i1,i2,i3,iv) + grad(i1,i2,i3,iv,3)*0.5
    rght(i1,i2,i3,iv) = prim(i1,i2,i3,iv) - grad(i1,i2,i3,iv,3)*0.5
  end do; end do; end do; end do
  !$omp end distribute parallel do simd
  !                  call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  !solverflux => flux(:,:,:,:,3)
  call solve (3, l+1, u, left, rght, flux,g1,g2)
  !                  call traceend (itimer1)
  !-----------------------------------------------------------------------------
  !call traceend (itimer)
END SUBROUTINE riemn3d

!===============================================================================
!> Update conservative variables with self%flux divergences, where self%flux(i1,:,:,:,1)
!> is centered at i1-0.5, so needs to be defined in the interval [li,uo]
!===============================================================================
SUBROUTINE divflux (self, mem, flux, nv, new, ds, dtime, l, u)
  !!$omp declare target
  class(solver_t),target:: self
  !.............................................................................
  integer:: i1, i2, i3, iv, dir, l(3), u(3)
  real:: dtd1, dtd2, dtd3
  integer, save:: itimer=0
  integer :: nv, new
  real, dimension(:,:,:,:,:,:),  pointer:: mem
  real, dimension(:,:,:,:,:), pointer :: flux
  real(8) :: dtime, ds(3)
  !mem => self%mem
  !flux => self%flux
  !nv = self%nv
  !new = self%new
  !-----------------------------------------------------------------------------
  !call tracebegin ('solver_t%divself%flux', itimer=itimer)
  dtd1 = dtime/ds(1)
  dtd2 = dtime/ds(2)
  dtd3 = dtime/ds(3)
  !l = self%li
  !u = self%ui
  !$omp distribute parallel do simd collapse(4)
  do iv=1,nv
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    mem(i1,i2,i3,iv,new,1) = &
    mem(i1,i2,i3,iv,new,1)   &
    - dtd1*(flux(i1+1,i2  ,i3  ,iv,1) - flux(i1,i2,i3,iv,1)) &
    - dtd2*(flux(i1  ,i2+1,i3  ,iv,2) - flux(i1,i2,i3,iv,2)) &
    - dtd3*(flux(i1  ,i2  ,i3+1,iv,3) - flux(i1,i2,i3,iv,3))
  end do
  end do
  end do
  end do
  !$omp end distribute parallel do simd
  !call traceend (itimer)
END SUBROUTINE divflux

!===============================================================================
!> Add source terms
!===============================================================================
SUBROUTINE sources (self, mem, prim, new, gamma)
  !!$omp declare target
  class(solver_t), target:: self
  !.............................................................................
  integer:: i1, i2, i3, iv, dir, new
  real:: dtd1, dtd2, dtd3
  real(8) :: gamma
  real, dimension(:,:,:,:), pointer :: prim
  real, dimension(:,:,:,:,:,:),  pointer:: mem
  real, dimension(:,:,:), pointer, contiguous:: p1, p2, p3, d, e
  logical:: forcing
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('solver_t%sources', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! If adjusting total energy, save kinetic part in self%prim(:,:,:,5)
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,1,new,1)
  p1 => mem(:,:,:,2,new,1)
  p2 => mem(:,:,:,3,new,1)
  p3 => mem(:,:,:,4,new,1)
  e  => mem(:,:,:,5,new,1)
  forcing = .false.
  !forcing = allocated(self%force_per_unit_mass) .or. &
  !!          allocated(self%force_per_unit_volume)
  if (gamma /= 1d0 .and. forcing) then
  !$omp parallel workshare
    prim(:,:,:,5) = p1**2 + p2**2 + p3**2
  !$omp end parallel workshare
  end if
  !-----------------------------------------------------------------------------
  ! Add forcing to momentum
  !-----------------------------------------------------------------------------
  !if (allocated(self%force_per_unit_mass)) then
  !  p1 = p1 + self%force_per_unit_mass(:,:,:,1)*d*self%dtime
  !  p2 = p2 + self%force_per_unit_mass(:,:,:,2)*d*self%dtime
  !  p3 = p3 + self%force_per_unit_mass(:,:,:,3)*d*self%dtime
  !end if
  !if (allocated(self%force_per_unit_volume)) then
  !  p1 = p1 + self%force_per_unit_volume(:,:,:,1)*self%dtime
  !  p2 = p2 + self%force_per_unit_volume(:,:,:,2)*self%dtime
  !  p3 = p3 + self%force_per_unit_volume(:,:,:,3)*self%dtime
  !end if
  !-----------------------------------------------------------------------------
  ! If adjusting total energy, use saved kinetic part
  !-----------------------------------------------------------------------------
  if (gamma == 1d0) then
  !$omp parallel workshare
    e = d*csound**2 + 0.5*(p1**2 + p2**2 + p3**2)/d
  !$omp end parallel workshare
  else if (forcing) then
  !$omp parallel workshare
    e = e + 0.5*(p1**2 + p2**2 + p3**2 - prim(:,:,:,5))/d
  !$omp end parallel workshare
  end if
  !call traceend (itimer)
END SUBROUTINE sources

!===============================================================================
!> Control timer verbosity; set detailed_timer=f to see only solver_t%update
!===============================================================================
SUBROUTINE tracebegin (label, itimer)
 character(len=*):: label
 integer:: itimer
 if (detailed_timer) then
 else
 end if
END SUBROUTINE tracebegin
!===============================================================================
SUBROUTINE traceend (itimer)
 integer:: itimer
 if (detailed_timer) then
 else
 end if
END SUBROUTINE traceend

END MODULE solver_mod
