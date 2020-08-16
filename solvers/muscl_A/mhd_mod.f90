!===============================================================================
!> MUSCL method using Riemann solver to compute fluxes
!>
!> The Riemann solver takes pointers to the left and right interface values
!> as input, derived from the primitive variables and their slopes.
!>
!> As primitive variables we choose log density, entropy per unit mass,
!> velocity and magnetic field. This ensures that the densities are positive, 
!> and that the slopes are derived from well-conditioned per-unit-mass variables.
!>
!> The number of guard zones needed is determined by, in the end, taking the
!> curl of E, which has a contribution eta*J, so J needs to have one guaerd zone.
!> Therefor, B has to have two guard zones, and A needs to have three.
!>
!> Variable no      conservative            primitive
!>      1               d                      ln(d)
!>      2               px                     ux
!>      3               pu                     ux
!>      4               pz                     ux
!>      5               s                      ss
!>      6               ax                     bx
!>      7               ay                     by
!>      8               az                     bz
!>
!>  Scratch:            name                   chunks
!>                      jj                     (:,:,:,3)
!>                      ee                     (:,:,:,3)
!>                      b_t                    (:,:,:,3)
!>                      grad                   (:,:,:,3,nv)
!>                      flux                   (:,:,:,nv,3)
!>                      ag                     (:,:,:,3,3)
!===============================================================================
MODULE mhd_mod
  USE io_mod
  USE trace_mod
  USE omp_lock_mod
  USE riemann_mod
  USE extras_mod
  USE scratch_mod
  implicit none
  private
  integer, parameter:: id=1, iu=2, iv=3, iw=4, is=5, ip=6
  type, public, extends(extras_t):: mhd_t
    real:: slope_type
    real:: eta
    real(8):: g1, g2
    type(riemann_t):: riemann
    logical:: mhd=.true.
  contains
    procedure:: init
    procedure:: update
  end type
  !-----------------------------------------------------------------------------
  ! Allocating & deallocating scratch arrays as part if hd_t carries a cost,
  ! (with gfortran the deallocate is more costly than the allocate), so a
  ! better alternative is to have permanent, thread-private scratch arrays:
  !-----------------------------------------------------------------------------
  real, dimension(:,:,:,:)  , allocatable, target:: prim            ! ,nv)
  real, dimension(:,:,:,:)  , allocatable        :: left            ! ,nv)
  real, dimension(:,:,:,:)  , allocatable        :: rght            ! ,nv)
  real, dimension(:,:,:,:)  , allocatable        :: jj              ! ,3)
  real, dimension(:,:,:,:)  , allocatable, target:: em              ! ,3)
  real, dimension(:,:,:,:)  , allocatable        :: b_t             ! ,3)
  real, dimension(:,:,:,:,:), allocatable        :: ag              ! ,3,3)
  real, dimension(:,:,:,:,:), allocatable, target:: grad            ! ,3,nv)
  real, dimension(:,:,:,:,:), allocatable        :: flux            ! ,nv,3)
  !$omp threadprivate (prim,left,rght,jj,em,b_t,ag,grad,flux)
  !-----------------------------------------------------------------------------
  ! Input parameters:
  !-----------------------------------------------------------------------------
  real:: gamma=1.4, csound=1.0, courant=0.2
  integer:: slope_type=2
  logical:: unsigned=.true., detailed_timer=.true.
  character(len=4):: method='hll'
  integer:: verbose=0
CONTAINS

!===============================================================================
!> Initialize parameters
!===============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  !.............................................................................
  namelist /muscl_params/ verbose, courant, gamma, csound, slope_type, &
    method, detailed_timer
  integer:: iostat, i
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%init')
  !call trace%print_id (id)
  !-----------------------------------------------------------------------------
  ! Ramses-like specific things
  !-----------------------------------------------------------------------------
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, muscl_params, iostat=iostat)
    call io%namelist_warning ('muscl_params', iostat)
    write (io%output, muscl_params)
    write (io%namelist, muscl_params)
  end if
  !$omp end critical (input_cr)
  self%riemann%method = method
  self%courant = courant
  self%csound = csound
  !-----------------------------------------------------------------------------
  ! Method depedent attributes
  !-----------------------------------------------------------------------------
  self%nw = 1
  if (trim(method) == 'hlld') then
    self%ng = 3
    !self%mhd = .true.
    self%kind = 'muscl_A_mhd_patch'
    if (self%nv == 0) self%nv = 9
  else
    self%ng = 2
    !self%mhd = .false.
    self%kind = 'muscl_A_hd_patch'
    if (self%nv == 0) self%nv = 6
  end if
  !-----------------------------------------------------------------------------
  ! Here we call patch_t%init directly, to initialize %mem and %mesh
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    write (io%output,'(a,2i4,l4)') 'mhd_t%init: nv, ng, mhd =', &
      self%nv, self%ng !, self%mhd
  end if
  call self%idx%init (mhd=.true., ie=5)
  call self%patch_t%init
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
  self%riemann%gamma = gamma
  self%riemann%g1 = self%g1
  self%riemann%g2 = self%g2
  !-----------------------------------------------------------------------------
  ! Generic settings
  !-----------------------------------------------------------------------------
  self%staggered = .false.
  do i=1,3
    self%mesh(i)%h = 0d0
  end do
  self%unsigned(self%idx%d) = unsigned
  self%unsigned(self%idx%s) = .false.
  self%pervolume(self%idx%px) = unsigned
  self%pervolume(self%idx%py) = unsigned
  self%pervolume(self%idx%pz) = unsigned
  !-----------------------------------------------------------------------------
  ! Calling gpatch_t%init, which may call ICs, must wait until all properties
  ! of patches (mesh, idx, patch_t%h staggering) are in place.
  !-----------------------------------------------------------------------------
  call self%gpatch_t%init
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Allocate work arrays
!===============================================================================
SUBROUTINE update (self)
  USE omp_timer_mod
  class(mhd_t):: self
  type(lock_t), pointer:: lock
  integer:: ng
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('mhd_t%update', itimer=itimer)
  call self%output
  call alloc   (self)
  lock => self%mem_lock(self%new)
  call ctoprim (self)
  if (verbose > 0) call self%stats (prim, 'prim')
  if (all(self%ng == 3)) then
    call slopes  (self, 2)
  else
    call slopes  (self, 1)
  end if
  if (verbose > 0) call self%stats (grad(:,:,:,:,1), 'grad')
  call predict (self)
  call riemn3d (self)
  if (verbose > 0) call self%stats (left, 'left')
  if (verbose > 0) call self%stats (rght, 'rght')
  if (verbose > 0) call self%stats (flux(:,:,:,:,1), 'flux')
  if (verbose > 0) print'("flux1:",1p,10e11.2)', flux(:,3,3,2,1)
  if (verbose > 0) print'("flux2:",1p,10e11.2)', flux(:,3,3,3,1)
  if (verbose > 0) print'("flux3:",1p,10e11.2)', flux(:,3,3,4,1)
  call divflux (self)
  if (verbose > 0) call self%stats (self%mem(:,:,:,:,self%new,1), 'mem')
  call sources (self)
  if (verbose > 0) call self%stats (self%mem(:,:,:,:,self%new,1), 'mem')
  call lock%unset ('update')
  call self%counter_update
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Allocate work arrays
!===============================================================================
SUBROUTINE alloc (self)
  class(mhd_t):: self
  integer, save:: itimer=0
  !.............................................................................
  call tracebegin ('mhd_t%alloc', itimer=itimer)
  if (.not.allocated(prim)) then
    allocate (prim (self%gn(1), self%gn(2), self%gn(3), self%nv   ))
    allocate (left (self%gn(1), self%gn(2) ,self%gn(3), self%nv   ))
    allocate (rght (self%gn(1), self%gn(2) ,self%gn(3), self%nv   ))
    if (self%mhd) then
    allocate (jj   (self%gn(1), self%gn(2) ,self%gn(3), 3         ))
    allocate (em   (self%gn(1), self%gn(2) ,self%gn(3), 3         ))
    allocate (b_t  (self%gn(1), self%gn(2) ,self%gn(3), 3         ))
    end if
    allocate (ag   (self%gn(1), self%gn(2) ,self%gn(3), 3      , 3))
    allocate (grad (self%gn(1), self%gn(2) ,self%gn(3), 3, self%nv))
    allocate (flux (self%gn(1), self%gn(2) ,self%gn(3), self%nv, 3))
    !---------------------------------------------------------------------------
    ! Set to zero once, to allow undefined guard zone values to be used
    !---------------------------------------------------------------------------
    prim = 0.0
    left = 0.0
    rght = 0.0
    grad = 0.0
    flux = 0.0
  end if
  call traceend (itimer)
END SUBROUTINE alloc

!===============================================================================
!> Deallocate work arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(mhd_t):: self
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('mhd_t%dealloc', itimer=itimer)
  deallocate (prim)
  deallocate (left)
  deallocate (rght)
  deallocate (grad)
  deallocate (flux)
  call trace%end (itimer)
END SUBROUTINE dealloc

!===============================================================================
!> Conservative to primitive transformation, using log(d), U, and Eint
!> Valid in full interval: [lb,ub] = [li-2,ui+2]
!===============================================================================
SUBROUTINE ctoprim (self)
  class(mhd_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, l(3), u(3), id, ip, is
  type(lock_t), pointer:: lock
  real:: u2, u2_max, c2, g2, b2, ss1
  real, contiguous, pointer, dimension(:,:,:)  :: d, s, ss
  real, contiguous, pointer, dimension(:,:,:,:):: a, b, uu
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('mhd_t%ctoprim', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Copy everything to new time slot, lock protected
  !-----------------------------------------------------------------------------
  lock => self%mem_lock(self%it)
  !call lock%set ('ctoprim')
  self%mem(:,:,:,:,self%new,1) = self%mem(:,:,:,:,self%it,1)
  !call lock%unset ('ctoprim')
  !-----------------------------------------------------------------------------
  ! Compute log density and velocities
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%new,1)
  s => self%mem(:,:,:,self%idx%s,self%new,1)
  ss => prim(:,:,:,self%idx%s)
  uu => prim(:,:,:,2:4)
  l = self%mesh%lb
  u = self%mesh%ub
  id = self%idx%d
  is = self%idx%s
  ip = self%idx%p
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    prim(i1,i2,i3,id) = log(self%mem(i1,i2,i3, 1,self%new,1)) ! log density
    ss1 = s(i1,i2,i3)/d(i1,i2,i3)
    prim(i1,i2,i3,is) = ss1                                   ! per unit mass
    prim(i1,i2,i3,ip) = &
      gamma*d(i1,i2,i3)*exp(ss(i1,i2,i3)*self%g1)             ! pressure
    prim(i1,i2,i3,2) = self%mem(i1,i2,i3,2,self%new,1) &      ! x-velocity
                     / self%mem(i1,i2,i3,1,self%new,1)
    prim(i1,i2,i3,3) = self%mem(i1,i2,i3,3,self%new,1) &      ! y-velocity
                     / self%mem(i1,i2,i3,1,self%new,1)
    prim(i1,i2,i3,4) = self%mem(i1,i2,i3,4,self%new,1) &      ! z-velocity
                     / self%mem(i1,i2,i3,1,self%new,1)
    prim(i1,i2,i3,5) = self%mem(i1,i2,i3,5,self%new,1) &      ! entropy / mass
                     / self%mem(i1,i2,i3,1,self%new,1)
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Compute B, J, and E
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    a => self%mem(:,:,:,self%idx%ax:self%idx%az,self%new,1)
    b => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%new,1)
    call curl (self, a, b, 2)
    call laplace3 (self, b, jj, 1)
    call cross (self, b, uu, em, 1)
    em = em + self%eta*jj
    call curl (self, em, b_t, 0)
  end if  
  !-----------------------------------------------------------------------------
  ! Compute self%u_max
  !-----------------------------------------------------------------------------
  u2_max = 0.0
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    if (self%mhd) then
      do i1=l(1),u(1)
        u2 = prim(i1,i2,i3,self%idx%px)**2 &
           + prim(i1,i2,i3,self%idx%py)**2 &
           + prim(i1,i2,i3,self%idx%pz)**2
        b2 = prim(i1,i2,i3,self%idx%bx)**2 &
           + prim(i1,i2,i3,self%idx%by)**2 &
           + prim(i1,i2,i3,self%idx%bz)**2
        c2 = (b2+gamma*prim(i1,i2,i3,ip))/d(i1,i2,i3)
        u2_max = max (u2_max, u2+c2)
      end do
    else
      do i1=l(1),u(1)
        u2 = prim(i1,i2,i3,self%idx%px)**2 &
           + prim(i1,i2,i3,self%idx%py)**2 &
           + prim(i1,i2,i3,self%idx%pz)**2
        c2 = gamma*prim(i1,i2,i3,ip)/d(i1,i2,i3)
        u2_max = max (u2_max, u2+c2)
      end do
    end if
  end do
  end do
  self%u_max = sqrt(u2_max)
  call self%patch_t%courant_condition
  call traceend (itimer)
END SUBROUTINE ctoprim

!===============================================================================
SUBROUTINE slopes (self, ng)
  class(mhd_t):: self
  integer:: ng
  !.............................................................................
  integer:: iv, l(3), u(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('mhd_t%slopes', itimer=itimer)
  l = self%mesh%li-ng
  u = self%mesh%ui+ng
!print '("slopes:",4(2x,3i3))', self%ng, l, u, ng
  if (slope_type == 1) then
    do iv=1,self%nv
      call minmod (self, prim(:,:,:,iv), grad(:,:,:,:,iv), l, u)
    end do
  else if (slope_type >= 2) then
    do iv=1,self%nv
      call moncen (self, prim(:,:,:,iv), grad(:,:,:,:,iv), l, u)
    end do
  else
    do iv=1,self%nv
      call central (self, prim(:,:,:,iv), grad(:,:,:,:,iv), l, u)
    end do
  end if
  call traceend (itimer)
END SUBROUTINE slopes

!===============================================================================
SUBROUTINE slope3 (self, pr, gr, ng)
  class(mhd_t):: self
  real, contiguous, dimension(:,:,:,:):: pr
  real, contiguous, dimension(:,:,:,:,:):: gr
  integer:: ng
  !.............................................................................
  integer:: i1, i2, i3, iv, dir, l(3), u(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('mhd_t%slope3', itimer=itimer)
  l = self%mesh%li-ng
  u = self%mesh%ui+ng
!print '("slope3:",4(2x,3i3))', self%ng, l, u, ng
  if (slope_type == 1) then
    do iv=1,3
      call minmod (self, pr(:,:,:,iv), gr(:,:,:,:,iv), l, u)
    end do
  else if (slope_type >= 2) then
    do iv=1,3
      call moncen (self, pr(:,:,:,iv), gr(:,:,:,:,iv), l, u)
    end do
  else
    do iv=1,3
      call central (self, pr(:,:,:,iv), gr(:,:,:,:,iv), l, u)
    end do
  end if
  call traceend (itimer)
END SUBROUTINE slope3

!===============================================================================
SUBROUTINE slope1 (self, pr, gr, ng)
  class(mhd_t):: self
  real, contiguous, dimension(:,:,:):: pr
  real, contiguous, dimension(:,:,:,:):: gr
  integer:: ng
  !.............................................................................
  integer:: i1, i2, i3, iv, dir, l(3), u(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('mhd_t%slope1', itimer=itimer)
  l = self%mesh%li-ng
  u = self%mesh%ui+ng
!print '("slope1:",4(2x,3i3))', self%ng, l, u, ng
  if (slope_type == 1) then
    call minmod (self, pr, gr, l, u)
  else if (slope_type >= 2) then
    call moncen (self, pr, gr, l, u)
  else
    call central (self, pr, gr, l, u)
  end if
  call traceend (itimer)
END SUBROUTINE slope1

!===============================================================================
!> Slope limiters, produce slopes in the interval [lb+1,ub-1] = [lo,uo]
!>
!> Costs at 16**3 point, in ms/call (divide with 4 to get mus/cell):
!>            T430
!> central    0.08
!> minmod     0.24
!> moncen     0.37
!===============================================================================
SUBROUTINE central (self, pr, gr, l, u)
  class(mhd_t):: self
  real, contiguous, dimension(:,:,:)  :: pr
  real, contiguous, dimension(:,:,:,:):: gr
  integer:: i1, i2, i3, l(3), u(3)
  !-----------------------------------------------------------------------------
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    do i1=l(1),u(1)
      gr (i1,i2,i3,1) = 0.5*(pr(i1+1,i2,i3) - pr(i1-1,i2,i3))
      gr (i1,i2,i3,2) = 0.5*(pr(i1,i2+1,i3) - pr(i1,i2-1,i3))
      gr (i1,i2,i3,3) = 0.5*(pr(i1,i2,i3+1) - pr(i1,i2,i3-1))
    end do
  end do
  end do
END SUBROUTINE central
!===============================================================================
SUBROUTINE minmod (self, pr, gr, l, u)
  class(mhd_t):: self
  real, contiguous, dimension(:,:,:)  :: pr
  real, contiguous, dimension(:,:,:,:):: gr
  real:: dlft, drgt
  integer:: i1, i2, i3, l(3), u(3)
  !-----------------------------------------------------------------------------
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !$omp simd
    do i1=l(1),u(1)
      ! slopes in first coordinate direction
      dlft = pr(i1  ,i2,i3) - pr(i1-1,i2,i3)
      drgt = pr(i1+1,i2,i3) - pr(i1  ,i2,i3)
      gr(i1,i2,i3,1) = merge(0.0,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=0.0)
      ! slopes in second coordinate direction
      dlft = pr(i1,i2  ,i3) - pr(i1,i2-1,i3)
      drgt = pr(i1,i2+1,i3) - pr(i1,i2  ,i3)
      gr(i1,i2,i3,2) = merge(0.0,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=0.0)
      ! slopes in third coordinate direction
      dlft = pr(i1,i2,i3  ) - pr(i1,i2,i3-1)
      drgt = pr(i1,i2,i3+1) - pr(i1,i2,i3  )
      gr(i1,i2,i3,3) = merge(0.0,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=0.0)
    end do
  end do
  end do
END SUBROUTINE minmod
!===============================================================================
SUBROUTINE moncen (self, pr, gr, l, u)
  class(mhd_t):: self
  real, contiguous, dimension(:,:,:)  :: pr
  real, contiguous, dimension(:,:,:,:):: gr
  real:: dlft,  drgt, dcen, dsgn, slop, dlim, c1, c2
  integer:: i1, i2, i3, l(3), u(3)
  !-----------------------------------------------------------------------------
  c2 = 0.5
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !$omp simd private(dlft,drgt,dcen,dsgn,slop,dlim)
    do i1 = l(1),u(1)
      ! slopes in first coordinate direction
      dlft = pr(i1  ,i2,i3) - pr(i1-1,i2,i3)
      drgt = pr(i1+1,i2,i3) - pr(i1  ,i2,i3)
      dcen = c2*(dlft+drgt)
      dsgn = sign(1.0, dcen)
      slop = min(abs(dlft),abs(drgt))
      dlim = merge(0.0,slop,(dlft*drgt)<=0.0)
      gr(i1,i2,i3,1) = dsgn*min(dlim,abs(dcen))
      ! slopes in second coordinate direction
      dlft = pr(i1,i2  ,i3) - pr(i1,i2-1,i3)
      drgt = pr(i1,i2+1,i3) - pr(i1,i2  ,i3)
      dcen = c2*(dlft+drgt)
      dsgn = sign(1.0,dcen)
      slop = min(abs(dlft),abs(drgt))
      dlim = merge(0.0,slop,(dlft*drgt)<=0.0)
      gr(i1,i2,i3,2) = dsgn*min(dlim,abs(dcen))
      ! slopes in third coordinate direction
      dlft = pr(i1,i2,i3  ) - pr(i1,i2,i3-1)
      drgt = pr(i1,i2,i3+1) - pr(i1,i2,i3  )
      dcen = c2*(dlft+drgt)
      dsgn = sign(1.0,dcen)
      slop = min(abs(dlft),abs(drgt))
      dlim = merge(0.0,slop,(dlft*drgt)<=0.0)
      gr(i1,i2,i3,3) = dsgn*min(dlim,abs(dcen))
    end do
  end do
  end do
END SUBROUTINE moncen

!===============================================================================
!> Predictor step for primitive variables; should produce left and right states
!> in intervals [li,uo] for left and [lo,uo] for right, so that after a shift
!> of right one index, the right state is valid in the interval [lo,ui]
!===============================================================================
SUBROUTINE predict (self)
  class(mhd_t):: self
  real, pointer, contiguous, dimension(:,:,:):: &
    ln_1, ln_2, ln_3, &
    ee_1, ee_2, ee_3, &
    ss_1, ss_2, ss_3, &
    u1_1, u1_2, u1_3, &
    u2_1, u2_2, u2_3, &
    u3_1, u3_2, u3_3, &
    ln,  ee, u1, u2, u3, d, ss, pp, pp_1, pp_2, pp_3
  real, pointer, contiguous, dimension(:,:,:,:):: uu, b
  real:: divu, dtd1, dtd2, dtd3, g1, dt2, lf1, lf2, lf3, qq, pod, em1, em2, em3
  real:: ff(size(prim,1),3)
  integer:: i1, i2, i3, l(3), u(3), ip
  integer, save:: itimer=0
  !.............................................................................
  call tracebegin ('mhd_t%predict', itimer=itimer)
  d => self%mem(:,:,:,self%idx%d,self%new,1)
  ip = self%idx%p
  if (self%nv >= 8) then
    b  => prim(:,:,:,6:8)
  end if
  uu   => prim(:,:,:,2:4)
  pp   => prim(:,:,:,ip)
  ln_1 => grad(:,:,:,1,1)
  ln_2 => grad(:,:,:,2,1)
  ln_3 => grad(:,:,:,3,1)
  u1_1 => grad(:,:,:,1,2)
  u1_2 => grad(:,:,:,2,2)
  u1_3 => grad(:,:,:,3,2)
  u2_1 => grad(:,:,:,1,3)
  u2_2 => grad(:,:,:,2,3)
  u2_3 => grad(:,:,:,3,3)
  u3_1 => grad(:,:,:,1,4)
  u3_2 => grad(:,:,:,2,4)
  u3_3 => grad(:,:,:,3,4)
  ss_1 => grad(:,:,:,1,5)
  ss_2 => grad(:,:,:,2,5)
  ss_3 => grad(:,:,:,3,5)
  pp_1 => grad(:,:,:,1,ip)
  pp_2 => grad(:,:,:,2,ip)
  pp_3 => grad(:,:,:,3,ip)
  ln   => prim(:,:,:,1)
  u1   => prim(:,:,:,2)
  u2   => prim(:,:,:,3)
  u3   => prim(:,:,:,4)
  ss   => prim(:,:,:,5)
  pp   => prim(:,:,:,ip)
  d    => self%mem(:,:,:,1,self%new,1)
  !-----------------------------------------------------------------------------
  ! Evolution of primitive variables over half a time step, with
  ! \partial\ln(d) / \partial t = -u*grad(ln_d) - div(u)
  !     \partial u / \partial t = -u*grad(u) - grad(P)/d
  !                     grad(P) = d*grad(e) + e*grad(d)
  !                             = d*grad(e) + e*d*grad(\ln(d))
  !     \partial S / \partial t = -u*grad(S) + Q_mass
  !     \partial B / \partial t = -curl(E)
  !                           E = eta*J - u x B
  !-----------------------------------------------------------------------------
  dt2  = self%dtime*0.5
  dtd1 = dt2/self%ds(1)
  dtd2 = dt2/self%ds(2)
  dtd3 = dt2/self%ds(3)
  g1 = self%g1
  l = self%mesh%lo
  u = self%mesh%uo
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    ff(:,:) = 0.0
    if (associated(scratch%force_per_unit_mass)) then
      ff(:,1) = ff(:,1) + scratch%force_per_unit_mass(:,i2,i3,1)
      ff(:,2) = ff(:,2) + scratch%force_per_unit_mass(:,i2,i3,2)
      ff(:,3) = ff(:,3) + scratch%force_per_unit_mass(:,i2,i3,3)
    end if
    if (associated(scratch%force_per_unit_volume)) then
      ff(:,1) = ff(:,1) + scratch%force_per_unit_volume(:,i2,i3,1)/d(:,i2,i3)
      ff(:,2) = ff(:,2) + scratch%force_per_unit_volume(:,i2,i3,2)/d(:,i2,i3)
      ff(:,3) = ff(:,3) + scratch%force_per_unit_volume(:,i2,i3,3)/d(:,i2,i3)
    end if
    if (self%mhd) then
      do i1=l(1),u(1)
        lf1 = jj(i1,i2,i3,2)*b(i1,i2,i3,3) - jj(i1,i2,i3,3)*b(i1,i2,i3,2)
        lf2 = jj(i1,i2,i3,3)*b(i1,i2,i3,1) - jj(i1,i2,i3,1)*b(i1,i2,i3,3)
        lf3 = jj(i1,i2,i3,1)*b(i1,i2,i3,2) - jj(i1,i2,i3,2)*b(i1,i2,i3,1)
        ff(i1,1) = ff(i1,1) + lf1
        ff(i1,2) = ff(i1,2) + lf2
        ff(i1,3) = ff(i1,3) + lf3
        em1 = b(i1,i2,i3,2)*uu(i1,i2,i3,3) - b(i1,i2,i3,3)*uu(i1,i2,i3,2)
        em2 = b(i1,i2,i3,3)*uu(i1,i2,i3,1) - b(i1,i2,i3,1)*uu(i1,i2,i3,3)
        em3 = b(i1,i2,i3,1)*uu(i1,i2,i3,2) - b(i1,i2,i3,2)*uu(i1,i2,i3,1)
        em(i1,i2,i3,1) = em1 + self%eta*jj(i2,i2,i3,3)
        em(i1,i2,i3,2) = em2 + self%eta*jj(i2,i2,i3,3)
        em(i1,i2,i3,3) = em3 + self%eta*jj(i2,i2,i3,3)
      end do
    end if
    do i1=l(1),u(1)
      divu = u1_1(i1,i2,i3)*dtd1 + u2_2(i1,i2,i3)*dtd2 + u3_3(i1,i2,i3)*dtd3
      ln(i1,i2,i3) = ln(i1,i2,i3) - ( &
        u1(i1,i2,i3)*ln_1(i1,i2,i3)*dtd1 + &
        u2(i1,i2,i3)*ln_2(i1,i2,i3)*dtd2 + &
        u3(i1,i2,i3)*ln_3(i1,i2,i3)*dtd3) - divu
      pod = prim(i1,i2,i3,ip)/d(i1,i2,i3)
      u1(i1,i2,i3) = u1(i1,i2,i3) - ( &
        u1(i1,i2,i3)*u1_1(i1,i2,i3)*dtd1 + &
        u2(i1,i2,i3)*u1_2(i1,i2,i3)*dtd2 + &
        u3(i1,i2,i3)*u1_3(i1,i2,i3)*dtd3) - &
        (ss_1(i1,i2,i3)*g1 + ln_1(i1,i2,i3)*gamma)*pod*dtd1 + &
        ff(i1,1)*dt2
      u2(i1,i2,i3) = u2(i1,i2,i3) - ( &
        u1(i1,i2,i3)*u2_1(i1,i2,i3)*dtd1 + &
        u2(i1,i2,i3)*u2_2(i1,i2,i3)*dtd2 + &
        u3(i1,i2,i3)*u2_3(i1,i2,i3)*dtd3) - &
        (ss_2(i1,i2,i3)*g1 + ln_2(i1,i2,i3)*gamma)*pod*dtd2 + &
        ff(i1,2)*dt2
      u3(i1,i2,i3) = u3(i1,i2,i3) - ( &
        u1(i1,i2,i3)*u3_1(i1,i2,i3)*dtd1 + &
        u2(i1,i2,i3)*u3_2(i1,i2,i3)*dtd2 + &
        u3(i1,i2,i3)*u3_3(i1,i2,i3)*dtd3) - &
        (ss_3(i1,i2,i3)*g1 + ln_3(i1,i2,i3)*gamma)*pod*dtd3 + &
        ff(i1,3)*dt2
      qq = 0.0
      ss(i1,i2,i3) = ss(i1,i2,i3) - ( &
        u1(i1,i2,i3)*ss_1(i1,i2,i3)*dtd1 + &
        u2(i1,i2,i3)*ss_2(i1,i2,i3)*dtd2 + &
        u3(i1,i2,i3)*ss_3(i1,i2,i3)*dtd3) + qq*dt2
    end do
  end do
  end do
  if (self%mhd) then
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      !-------------------------------------------------------------------------
      ! \partial B/\partial t = -curl(E)
      !-------------------------------------------------------------------------
      b(i1,i2,i3,1) = b(i1,i2,i3,1) + b_t(i1,i2,i3,1)*dt2
      b(i1,i2,i3,1) = b(i1,i2,i3,1) + b_t(i1,i2,i3,1)*dt2
      b(i1,i2,i3,1) = b(i1,i2,i3,1) + b_t(i1,i2,i3,1)*dt2
    end do
    end do
    end do
  end if
  call traceend (itimer)
END SUBROUTINE predict

!===============================================================================
!> Select the left and right states and call the Riemann solver to get fluxes.
!> The left and right states are defined for indices [lo,uo], which is where
!> the slopes are valid.  The left states refer to positions [lo+0.5,uo+0.5]
!> and the right states refer to positions [lo-0.5,uo-0.5], while we need the
!> fluxes in the interval [li-0.5,ui+0.5] = [lo+0.5,uo-0.5], which matches
!> what the left and right states combined can supply.  In order to refer to
!> left and right states at the same physical location, we stay with the
!> indices of the right states, whose physical locations match those where we
!> want the fluxes, while referring to the left state one index to the left
!> inside the Riemann solver
!===============================================================================
SUBROUTINE riemn3d (self)
  class(mhd_t):: self
  !.............................................................................
  integer:: l(3), u(3), i1, i2, i3, iv
  integer, save:: itimer=0, itimer1=0
  !-----------------------------------------------------------------------------
  call tracebegin ('mhd_t%riemn3d', itimer=itimer)
  l = self%mesh%lo
  u = self%mesh%uo
  !-----------------------------------------------------------------------------
  ! X-direction
  !-----------------------------------------------------------------------------
  do iv=1,self%nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    left(i1,i2,i3,iv) = prim(i1,i2,i3,iv) + grad(i1,i2,i3,1,iv)*0.5
    rght(i1,i2,i3,iv) = prim(i1,i2,i3,iv) - grad(i1,i2,i3,1,iv)*0.5
  end do; end do; end do; end do
                    call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  call self%riemann%solve (1, l+1, u, left, rght, flux(:,:,:,:,1))
                    call traceend (itimer1)
  !-----------------------------------------------------------------------------
  ! Y-direction
  !-----------------------------------------------------------------------------
  do iv=1,self%nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    left(i1,i2,i3,iv) = prim(i1,i2,i3,iv) + grad(i1,i2,i3,2,iv)*0.5
    rght(i1,i2,i3,iv) = prim(i1,i2,i3,iv) - grad(i1,i2,i3,2,iv)*0.5
  end do; end do; end do; end do
                    call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  call self%riemann%solve (2, l+1, u, left, rght, flux(:,:,:,:,2))
                    call traceend (itimer1)
  !-----------------------------------------------------------------------------
  ! Z-direction
  !-----------------------------------------------------------------------------
  do iv=1,self%nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    left(i1,i2,i3,iv) = prim(i1,i2,i3,iv) + grad(i1,i2,i3,3,iv)*0.5
    rght(i1,i2,i3,iv) = prim(i1,i2,i3,iv) - grad(i1,i2,i3,3,iv)*0.5
  end do; end do; end do; end do
                    call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  call self%riemann%solve (3, l+1, u, left, rght, flux(:,:,:,:,3))
                    call traceend (itimer1)
  !-----------------------------------------------------------------------------
  call traceend (itimer)
END SUBROUTINE riemn3d

!===============================================================================
!> Update conservative variables with flux divergences, where flux(i1,:,:,:,1)
!> is centered at i1-0.5, so needs to be defined in the interval [li,uo]
!===============================================================================
SUBROUTINE divflux (self)
  class(mhd_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, dir, l(3), u(3)
  real:: dtd1, dtd2, dtd3
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('mhd_t%divflux', itimer=itimer)
  dtd1 = self%dtime/self%ds(1)
  dtd2 = self%dtime/self%ds(2)
  dtd3 = self%dtime/self%ds(3)
  l = self%mesh%li
  u = self%mesh%ui
  do iv=1,self%nv
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    self%mem(i1,i2,i3,iv,self%new,1) = &
    self%mem(i1,i2,i3,iv,self%new,1)   &
    - dtd1*(flux(i1+1,i2  ,i3  ,iv,1) - flux(i1,i2,i3,iv,1)) &
    - dtd2*(flux(i1  ,i2+1,i3  ,iv,2) - flux(i1,i2,i3,iv,2)) &
    - dtd3*(flux(i1  ,i2  ,i3+1,iv,3) - flux(i1,i2,i3,iv,3))
  end do
  end do
  end do
  end do
  call traceend (itimer)
END SUBROUTINE divflux

!===============================================================================
!> Add source terms
!===============================================================================
SUBROUTINE sources (self)
  class(mhd_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, dir
  real:: dtd1, dtd2, dtd3
  real, dimension(:,:,:), pointer, contiguous:: p1, p2, p3, d, s, b1, b2, b3
  logical:: forcing
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('mhd_t%sources', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! If adjusting total energy, save kinetic part in prim(:,:,:,5)
  !-----------------------------------------------------------------------------
  d  => self%mem(:,:,:,1,self%new,1)
  p1 => self%mem(:,:,:,2,self%new,1)
  p2 => self%mem(:,:,:,3,self%new,1)
  p3 => self%mem(:,:,:,4,self%new,1)
  forcing = associated(scratch%force_per_unit_mass) .or. &
            associated(scratch%force_per_unit_volume)
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    b1 => self%mem(:,:,:,6,self%new,1)
    b2 => self%mem(:,:,:,7,self%new,1)
    b3 => self%mem(:,:,:,8,self%new,1)
  end if
  if (self%gamma /= 1d0 .and. forcing) then
    if (self%mhd) then
      b1 => self%mem(:,:,:,6,self%new,1)
      b2 => self%mem(:,:,:,7,self%new,1)
      b3 => self%mem(:,:,:,8,self%new,1)
      prim(:,:,:,5) = (p1**2 + p2**2 + p3**2)/d + (b1**2 + b2**2 + b3**2)
    else
      prim(:,:,:,5) = p1**2 + p2**2 + p3**2
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Add forcing to momentum
  !-----------------------------------------------------------------------------
  if (associated(scratch%force_per_unit_mass)) then
    p1 = p1 + scratch%force_per_unit_mass(:,:,:,1)*d*self%dtime
    p2 = p2 + scratch%force_per_unit_mass(:,:,:,2)*d*self%dtime
    p3 = p3 + scratch%force_per_unit_mass(:,:,:,3)*d*self%dtime
  end if
  if (associated(scratch%force_per_unit_volume)) then
    p1 = p1 + scratch%force_per_unit_volume(:,:,:,1)*self%dtime
    p2 = p2 + scratch%force_per_unit_volume(:,:,:,2)*self%dtime
    p3 = p3 + scratch%force_per_unit_volume(:,:,:,3)*self%dtime
  end if
  call traceend (itimer)
END SUBROUTINE sources

!===============================================================================
!> Central difference curl operator
!===============================================================================
SUBROUTINE cross (self, a, b, c, ng)
  class(mhd_t):: self
  real, contiguous, dimension(:,:,:,:):: a, b, c
  integer:: ng
  !.............................................................................
  integer:: i1, i2, i3
  !-----------------------------------------------------------------------------
  do i3=self%mesh(3)%li+ng,self%mesh(3)%ui+ng  
  do i2=self%mesh(2)%li+ng,self%mesh(2)%ui+ng
  do i1=self%mesh(1)%li+ng,self%mesh(1)%ui+ng
    c(i1,i2,i3,1) = a(i1,i2,i3,2)*b(i1,i2,i3,3) &
                  - a(i1,i2,i3,3)*b(i1,i2,i3,2)
    c(i1,i2,i3,2) = a(i1,i2,i3,3)*b(i1,i2,i3,1) &
                  - a(i1,i2,i3,1)*b(i1,i2,i3,3)
    c(i1,i2,i3,3) = a(i1,i2,i3,1)*b(i1,i2,i3,2) &
                  - a(i1,i2,i3,2)*b(i1,i2,i3,1)
  end do
  end do
  end do
END SUBROUTINE cross

!===============================================================================
!> Central difference curl operator
!===============================================================================
SUBROUTINE curl (self, a, b, ng)
  class(mhd_t):: self
  real, contiguous, dimension(:,:,:,:):: a, b
  integer:: ng
  !.............................................................................
  integer:: i1, i2, i3
  real:: c(3)
  !-----------------------------------------------------------------------------
  call slope3 (self, a, ag, ng)
  c = 1./self%ds
  do i3=self%mesh(3)%li+ng,self%mesh(3)%ui+ng  
  do i2=self%mesh(2)%li+ng,self%mesh(2)%ui+ng
  do i1=self%mesh(1)%li+ng,self%mesh(1)%ui+ng
    b(i1,i2,i3,1) = ag(i1,i2,i3,3,2)-ag(i1,i2,i3,2,3)
    b(i1,i2,i3,2) = ag(i1,i2,i3,1,3)-ag(i1,i2,i3,3,1)
    b(i1,i2,i3,3) = ag(i1,i2,i3,2,1)-ag(i1,i2,i3,1,2)
  end do
  end do
  end do
END SUBROUTINE curl

!===============================================================================
!> Central difference Laplace operator
!===============================================================================
SUBROUTINE laplace3 (self, a, b, ng)
  class(mhd_t):: self
  real, contiguous, dimension(:,:,:,:):: a, b
  integer:: ng
  !.............................................................................
  integer:: i1, i2, i3, dir
  real:: c(3)
  !-----------------------------------------------------------------------------
  c = 1./self%ds**2
  do dir=1,3
  do i3=self%mesh(3)%li-ng,self%mesh(3)%ui+ng
  do i2=self%mesh(2)%li-ng,self%mesh(2)%ui+ng
  do i1=self%mesh(1)%li-ng,self%mesh(1)%ui+ng
    b(i1,i2,i3,dir) = a(i1+1,i2  ,i3  ,dir)+a(i1-1,i2  ,i3  ,dir)*c(1) &
                    + a(i1  ,i2+1,i3  ,dir)+a(i1  ,i2-1,i3  ,dir)*c(2) &
                    + a(i1  ,i2  ,i3+1,dir)+a(i1  ,i2  ,i3-1,dir)*c(3) &
                    - a(i1  ,i2  ,i3  ,dir)*(c(1)+c(2)+c(3))*2.0
  end do
  end do
  end do
  end do
END SUBROUTINE laplace3

!===============================================================================
!> Control timer verbosity; set detailed_timer=f to see only hd_t%update
!===============================================================================
SUBROUTINE tracebegin (label, itimer)
 character(len=*):: label
 integer:: itimer
 if (detailed_timer) then
   call trace%begin (label, itimer=itimer)
 else
   call trace%begin (label)
 end if
END SUBROUTINE tracebegin
!===============================================================================
SUBROUTINE traceend (itimer)
 integer:: itimer
 if (detailed_timer) then
   call trace%end (itimer)
 else
   call trace%end ()
 end if
END SUBROUTINE traceend

END MODULE mhd_mod
