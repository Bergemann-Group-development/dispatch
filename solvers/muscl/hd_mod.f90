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
MODULE hd_mod
  USE io_mod
  USE trace_mod
  USE omp_lock_mod
  USE riemann_mod
  USE extras_mod
  USE scratch_mod
  implicit none
  private
  type, public, extends(extras_t):: hd_t
    real:: slope_type
    real(8):: g1, g2
    type(riemann_t):: riemann
  contains
    procedure:: init
    procedure:: update
  end type
  !-----------------------------------------------------------------------------
  ! Allocating & deallocating scratch arrays as part if hd_t carries a cost,
  ! (with gfortran the deallocate is more costly than the allocate), so a
  ! better alternative is to have permanent, thread-private scratch arrays:
  !-----------------------------------------------------------------------------
  real, dimension(:,:,:,:)  , pointer, contiguous:: prim => null()
  real, dimension(:,:,:,:)  , pointer, contiguous:: left => null()
  real, dimension(:,:,:,:)  , pointer, contiguous:: rght => null()
  real, dimension(:,:,:,:,:), pointer, contiguous:: grad => null()
  real, dimension(:,:,:,:,:), pointer, contiguous:: flux => null()
  !$omp threadprivate (prim,left,rght,grad,flux)
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
  class(hd_t):: self
  !.............................................................................
  namelist /muscl_params/ verbose, courant, gamma, csound, slope_type, &
    method, detailed_timer
  integer:: iostat, i
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('hd_t%init')
  !call trace%print_id (id)
  self%nw = 1
  self%ng = 2
  self%kind = 'ramses_hd_patch'
  if (self%nv == 0) self%nv = 5
  !-----------------------------------------------------------------------------
  ! Here we call patch_t%init directly, to initialize %mem and %mesh
  !-----------------------------------------------------------------------------
  call self%idx%init (mhd=.false., ienergy=5)
  call self%patch_t%init
  !-----------------------------------------------------------------------------
  ! Ramses-like specific things
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, muscl_params, iostat=iostat)
    call io%namelist_warning ('muscl_params', iostat)
    write (io%output, muscl_params)
    write (io%namelist, muscl_params)
  end if
  call io%lock%unset
  self%riemann%method = method
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
  self%unsigned(self%idx%e) = unsigned
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
  class(hd_t):: self
  integer, save:: itimer=0
  type(lock_t), pointer:: lock
  !.............................................................................
  call trace%begin ('hd_t%update', itimer=itimer)
  call self%output
  call alloc   (self)
  lock => self%mem_lock(self%new)
  call ctoprim (self)
  call slopes  (self)
  call predict (self)
  call riemn3d (self)
  call divflux (self)
  call sources (self)
  call lock%unset ('update')
  call self%counter_update
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Allocate work arrays
!===============================================================================
SUBROUTINE alloc (self)
  class(hd_t):: self
  integer, save:: itimer=0
  !.............................................................................
  call tracebegin ('hd_t%alloc', itimer=itimer)
  if (.not.associated(prim)) then
    allocate (prim (self%gn(1), self%gn(2), self%gn(3), self%nv   ))
    allocate (left (self%gn(1), self%gn(2) ,self%gn(3), self%nv   ))
    allocate (rght (self%gn(1), self%gn(2) ,self%gn(3), self%nv   ))
    allocate (grad (self%gn(1), self%gn(2) ,self%gn(3), self%nv, 3))
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
  class(hd_t):: self
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('hd_t%dealloc', itimer=itimer)
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
  class(hd_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, l(3), u(3)
  type(lock_t), pointer:: lock
  real:: u2, u2_max, c2, g2
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('hd_t%ctoprim', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Copy everything to new time slot, lock protected
  !-----------------------------------------------------------------------------
  lock => self%mem_lock(self%it)
  !call lock%set ('ctoprim')
  self%mem(:,:,:,:,self%new,1) = self%mem(:,:,:,:,self%it,1)
  !call lock%unset ('ctoprim')
  !-----------------------------------------------------------------------------
  l = self%mesh%lb
  u = self%mesh%ub
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    prim(:,i2,i3,1) = log(self%mem(:,i2,i3, 1,self%new,1))
    do iv=2,5
      prim(:,i2,i3,iv)  = self%mem(:,i2,i3,iv,self%new,1) &
                        / self%mem(:,i2,i3, 1,self%new,1)
    end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Subtract kinetic energy per unit mass and compute self%u_max
  !-----------------------------------------------------------------------------
  if (self%gamma == 1d0) then
    g2 = 1.0
  else
    g2 = (self%gamma-1d0)*self%gamma
  end if
  u2_max = 0.0
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
  self%u_max = sqrt(u2_max)
!print *,self%id,self%u_max
  call self%patch_t%courant_condition
  call traceend (itimer)
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
SUBROUTINE slopes (self)
  class(hd_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, dir, l(3), u(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('hd_t%slopes', itimer=itimer)
  l = self%mesh%lb+1
  u = self%mesh%ub-1
  if (slope_type == 1) then
    call minmod
  else if (slope_type >= 2) then
    call moncen
  else
    call central
  end if
  call traceend (itimer)
contains
!===============================================================================
subroutine central
  do iv=1,self%nv
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
end subroutine central
!===============================================================================
subroutine minmod
  real:: dlft, drgt
  do iv=1,self%nv
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !$omp simd
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
end subroutine minmod
!===============================================================================
subroutine moncen
  real:: dlft,  drgt, dcen, dsgn, slop, dlim, c1, c2
  c2 = 0.5
  do iv=1,self%nv
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
end subroutine moncen
END SUBROUTINE slopes

!===============================================================================
!> Predictor step for primitive variables; should produce left and right states
!> in intervals [li,uo] for left and [lo,uo] for right, so that after a shift
!> of right one index, the right state is valid in the interval [lo,ui]
!===============================================================================
SUBROUTINE predict (self)
  class(hd_t):: self
  real, pointer, contiguous, dimension(:,:,:):: &
    ln_1, ln_2, ln_3, &
    ee_1, ee_2, ee_3, &
    u1_1, u1_2, u1_3, &
    u2_1, u2_2, u2_3, &
    u3_1, u3_2, u3_3, &
    ln, ee, u1, u2, u3, d
  real:: divu, dtd1, dtd2, dtd3, g1, dt2
  real:: ff(size(prim,1),3)
  integer:: i1, i2, i3, l(3), u(3)
  integer, save:: itimer=0
  !.............................................................................
  call tracebegin ('hd_t%predict', itimer=itimer)
  ln_1 => grad(:,:,:,1,1)
  ln_2 => grad(:,:,:,1,2)
  ln_3 => grad(:,:,:,1,3)
  u1_1 => grad(:,:,:,2,1)
  u1_2 => grad(:,:,:,2,2)
  u1_3 => grad(:,:,:,2,3)
  u2_1 => grad(:,:,:,3,1)
  u2_2 => grad(:,:,:,3,2)
  u2_3 => grad(:,:,:,3,3)
  u3_1 => grad(:,:,:,4,1)
  u3_2 => grad(:,:,:,4,2)
  u3_3 => grad(:,:,:,4,3)
  ee_1 => grad(:,:,:,5,1)
  ee_2 => grad(:,:,:,5,2)
  ee_3 => grad(:,:,:,5,3)
  ln   => prim(:,:,:,1)
  u1   => prim(:,:,:,2)
  u2   => prim(:,:,:,3)
  u3   => prim(:,:,:,4)
  ee   => prim(:,:,:,5)
  d    => self%mem(:,:,:,1,self%new,1)
  !-----------------------------------------------------------------------------
  ! Evolution of primitive variables over half a time step, with
  ! \partial\ln(d) / \partial t = -u*grad(ln_d) - div(u)
  !     \partial u / \partial t = -u*grad(u) - grad(P)/d
  !                     grad(P) = d*grad(e) + e*grad(d)
  !                             = d*grad(e) + e*d*grad(\ln(d))
  !     \partial e / \partial t = -u*grad(e) - (P/d)*div(U)
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
    do i1=l(1),u(1)
      divu = u1_1(i1,i2,i3)*dtd1 + u2_2(i1,i2,i3)*dtd2 + u3_3(i1,i2,i3)*dtd3
      ln(i1,i2,i3) = ln(i1,i2,i3) - ( &
        u1(i1,i2,i3)*ln_1(i1,i2,i3)*dtd1 + &
        u2(i1,i2,i3)*ln_2(i1,i2,i3)*dtd2 + &
        u3(i1,i2,i3)*ln_3(i1,i2,i3)*dtd3) - divu
      u1(i1,i2,i3) = u1(i1,i2,i3) - ( &
        u1(i1,i2,i3)*u1_1(i1,i2,i3)*dtd1 + &
        u2(i1,i2,i3)*u1_2(i1,i2,i3)*dtd2 + &
        u3(i1,i2,i3)*u1_3(i1,i2,i3)*dtd3) - &
        (ee_1(i1,i2,i3) + ee(i1,i2,i3)*ln_1(i1,i2,i3))*g1*dtd1 + ff(i1,1)*dt2
      u2(i1,i2,i3) = u2(i1,i2,i3) - ( &
        u1(i1,i2,i3)*u2_1(i1,i2,i3)*dtd1 + &
        u2(i1,i2,i3)*u2_2(i1,i2,i3)*dtd2 + &
        u3(i1,i2,i3)*u2_3(i1,i2,i3)*dtd3) - &
        (ee_2(i1,i2,i3) + ee(i1,i2,i3)*ln_2(i1,i2,i3))*g1*dtd2 + ff(i1,2)*dt2
      u3(i1,i2,i3) = u3(i1,i2,i3) - ( &
        u1(i1,i2,i3)*u3_1(i1,i2,i3)*dtd1 + &
        u2(i1,i2,i3)*u3_2(i1,i2,i3)*dtd2 + &
        u3(i1,i2,i3)*u3_3(i1,i2,i3)*dtd3) - &
        (ee_3(i1,i2,i3) + ee(i1,i2,i3)*ln_3(i1,i2,i3))*g1*dtd3 + ff(i1,3)*dt2
      ee(i1,i2,i3) = ee(i1,i2,i3) - ( &
        u1(i1,i2,i3)*ee_1(i1,i2,i3)*dtd1 + &
        u2(i1,i2,i3)*ee_2(i1,i2,i3)*dtd2 + &
        u3(i1,i2,i3)*ee_3(i1,i2,i3)*dtd3) - divu*g1
    end do
  end do
  end do
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
  class(hd_t):: self
  !.............................................................................
  integer:: l(3), u(3), i1, i2, i3, iv
  integer, save:: itimer=0, itimer1=0
  !-----------------------------------------------------------------------------
  call tracebegin ('hd_t%riemn3d', itimer=itimer)
  l = self%mesh%lo
  u = self%mesh%uo
  !-----------------------------------------------------------------------------
  ! X-direction
  !-----------------------------------------------------------------------------
  do iv=1,self%nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    left(i1,i2,i3,iv) = prim(i1,i2,i3,iv) + grad(i1,i2,i3,iv,1)*0.5
    rght(i1,i2,i3,iv) = prim(i1,i2,i3,iv) - grad(i1,i2,i3,iv,1)*0.5
  end do; end do; end do; end do
                    call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  call self%riemann%solve (1, l+1, u, left, rght, flux(:,:,:,:,1))
                    call traceend (itimer1)
  !-----------------------------------------------------------------------------
  ! Y-direction
  !-----------------------------------------------------------------------------
  do iv=1,self%nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    left(i1,i2,i3,iv) = prim(i1,i2,i3,iv) + grad(i1,i2,i3,iv,2)*0.5
    rght(i1,i2,i3,iv) = prim(i1,i2,i3,iv) - grad(i1,i2,i3,iv,2)*0.5
  end do; end do; end do; end do
                    call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  call self%riemann%solve (2, l+1, u, left, rght, flux(:,:,:,:,2))
                    call traceend (itimer1)
  !-----------------------------------------------------------------------------
  ! Z-direction
  !-----------------------------------------------------------------------------
  do iv=1,self%nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    left(i1,i2,i3,iv) = prim(i1,i2,i3,iv) + grad(i1,i2,i3,iv,3)*0.5
    rght(i1,i2,i3,iv) = prim(i1,i2,i3,iv) - grad(i1,i2,i3,iv,3)*0.5
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
  class(hd_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, dir, l(3), u(3)
  real:: dtd1, dtd2, dtd3
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('hd_t%divflux', itimer=itimer)
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
  class(hd_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, dir
  real:: dtd1, dtd2, dtd3
  real, dimension(:,:,:), pointer, contiguous:: p1, p2, p3, d, e
  logical:: forcing
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call tracebegin ('hd_t%sources', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! If adjusting total energy, save kinetic part in prim(:,:,:,5)
  !-----------------------------------------------------------------------------
  d  => self%mem(:,:,:,1,self%new,1)
  p1 => self%mem(:,:,:,2,self%new,1)
  p2 => self%mem(:,:,:,3,self%new,1)
  p3 => self%mem(:,:,:,4,self%new,1)
  e  => self%mem(:,:,:,5,self%new,1)
  forcing = associated(scratch%force_per_unit_mass) .or. &
            associated(scratch%force_per_unit_volume)
  if (self%gamma /= 1d0 .and. forcing) then
    prim(:,:,:,5) = p1**2 + p2**2 + p3**2
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
  !-----------------------------------------------------------------------------
  ! If adjusting total energy, use saved kinetic part
  !-----------------------------------------------------------------------------
  if (self%gamma == 1d0) then
    e = d*csound**2 + 0.5*(p1**2 + p2**2 + p3**2)/d
  else if (forcing) then
    e = e + 0.5*(p1**2 + p2**2 + p3**2 - prim(:,:,:,5))/d
  end if
  call traceend (itimer)
END SUBROUTINE sources

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

END MODULE hd_mod
