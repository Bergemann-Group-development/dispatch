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
  !USE io_mod
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
    real, pointer:: mem(:,:,:,:,:,:)
    real, dimension(:,:,:,:)  , pointer, contiguous:: prim
    real, dimension(:,:,:,:)  , pointer, contiguous:: left
    real, dimension(:,:,:,:)  , pointer, contiguous:: rght
    real, dimension(:,:,:,:,:), pointer, contiguous:: grad
    real, dimension(:,:,:,:,:), pointer, contiguous:: flux
  contains
    procedure:: init
    procedure:: update
  end type
  real(8):: start
  !-----------------------------------------------------------------------------
  ! Allocating & deallocating scratch arrays as part if solver_t carries a cost,
  ! (with gfortran the deallocate is more costly than the allocate), so a
  ! better alternative is to have permanent, thread-private scratch arrays:
  !-----------------------------------------------------------------------------
  !!omp threadprivate (prim,left,rght,grad,flux)
  !-----------------------------------------------------------------------------
  ! Input parameters:
  !-----------------------------------------------------------------------------
  integer:: n=NSIZE                        ! internal cells
  integer:: ng=2                        ! guard cells
  integer:: slope_type=2                ! slope limiter
  integer:: verbose=0                   ! info level
  logical:: detailed_timer=.true.       ! optional timer
  real(8):: gamma=1.4                   ! gas gamma
  real:: csound=1.0                     ! sound speed
  real:: courant=0.2                    ! courant number
  character(len=4):: method='hll'       ! riemann method
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
  !integer:: m
  integer:: m, i1, i2, i3
  real(8):: r1, r2, r3
  real(8), parameter:: pi2=8d0*atan(1d0)
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
  !self%mem(:,:,:,1,:,:) = 1.0
  !self%mem(:,:,:,2,:,:) = 0.0
  !self%mem(:,:,:,3,:,:) = 0.0
  !self%mem(:,:,:,4,:,:) = 0.0
  !self%mem(:,:,:,5,:,:) = 1.0
  self%mem(:,:,:,1,:,:) = 1.0
  self%mem(:,:,:,5,:,:) = 1.0
  do i3=1,self%gn(3)
    r3 = (i3-1)/self%gn(3)
    do i2=1,self%gn(2)
      r2 = (i2-1)/self%gn(2)
      do i1=1,self%gn(1)
        r1 = (i1-1)/self%gn(1)
        self%mem(i1,i2,i3,2,:,:) = sin(r1*pi2)
        self%mem(i1,i2,i3,3,:,:) = sin(r2*pi2)
        self%mem(i1,i2,i3,4,:,:) = sin(r3*pi2)
      end do
    end do
  end do
  !call io%lock%set
  if (first_time) then
#ifdef DBG     
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
#endif    
    start = omp_get_wtime()
    first_time = .false.
  end if
  !call io%lock%unset
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

SUBROUTINE gpu_work (end_time, self_nv, self_gn, self_mem, self_it, self_new, local_slope_type, local_csound, self_time, self_dtime, self_lb, self_ub, self_gamma, self_ds, self_g1, self_lo, self_uo, self_g2, self_li, self_ui, self_prim, self_left, self_rght, self_grad, self_flux)
  implicit none
  real(8):: end_time
  integer:: self_nv
  integer, dimension(3):: self_gn, self_lb, self_ub, self_lo, self_uo, self_li, self_ui
  real, pointer:: self_mem(:,:,:,:,:,:)
  integer:: self_it, self_new
  integer:: local_slope_type 
  real:: local_csound
  real(8):: self_time, self_dtime, self_gamma, self_g1, self_g2
  real(8):: self_ds(3) 

  logical:: continue_running
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_left
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_rght
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_grad
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_flux
  !$omp declare target

  continue_running = .true.
  
  do while (continue_running)
    call ctoprim (self_it, self_new, self_lb, self_ub, self_mem, self_prim, self_gamma)
    call slopes  (self_nv, self_lb, self_ub, local_slope_type, self_prim, self_grad)
    call predict (self_g1, self_lo, self_uo, self_dtime, self_new, self_ds, self_prim, self_grad, self_mem)
    call riemn3d (self_g1, self_g2, self_lo, self_uo, self_nv, self_left, self_rght, self_prim, self_grad, self_flux)
    call divflux (self_new, self_dtime, self_ds, self_li, self_ui, self_nv, self_mem, self_flux)
    call sources (local_csound, self_new, self_gamma, self_mem, self_prim)
    self_time = self_time + self_dtime
    self_new = mod(self_new, 5)+1
    self_it  = mod(self_it,  5)+1
    continue_running = self_time .lt.  end_time
  end do

END SUBROUTINE gpu_work

SUBROUTINE openmp_work (end_time, self_nv, self_gn, self_mem, self_it, self_new, local_slope_type, local_csound, self_time, self_dtime, self_lb, self_ub, self_gamma, self_ds, self_g1, self_lo, self_uo, self_g2, self_li, self_ui)
  implicit none
  real(8):: end_time
  integer:: self_nv
  integer, dimension(3):: self_gn, self_lb, self_ub, self_lo, self_uo, self_li, self_ui
  real, pointer:: self_mem(:,:,:,:,:,:)
  integer:: self_it, self_new
  integer:: local_slope_type 
  real:: local_csound
  real(8):: self_time, self_dtime, self_gamma, self_g1, self_g2
  real(8):: self_ds(3) 

  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_left
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_rght
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_grad
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_flux
  
  call alloc   (self_gn, self_nv, self_prim, self_left, self_rght, self_grad, self_flux)

  !$omp target map(to: self_gn, self_nv, end_time, local_slope_type, local_csound, self_dtime, self_lo, self_uo, self_g1, self_g2, self_ds, self_li, self_ui, self_prim, self_left, self_rght, self_grad, self_flux) map(tofrom: self_mem, self_it, self_new, self_time) defaultmap(tofrom:scalar)

  call gpu_work (end_time, self_nv, self_gn, self_mem, self_it, self_new, local_slope_type, local_csound, self_time, self_dtime, self_lb, self_ub, self_gamma, self_ds, self_g1, self_lo, self_uo, self_g2, self_li, self_ui, self_prim, self_left, self_rght, self_grad, self_flux)
  
  !$omp end target

  call dealloc (self_prim, self_left, self_rght, self_grad, self_flux)

END SUBROUTINE openmp_work

!===============================================================================
!> Allocate work arrays
!===============================================================================
SUBROUTINE update (self, end_time)
  implicit none
  class(solver_t):: self
  real(8):: end_time

  ! integer, save:: itimer=0
  real(8):: time

  ! to be transferred to device
  integer:: self_it, self_new ! tofrom
  real(8):: self_time ! tofrom
  !.............................................................................
  !call self%output

  time = omp_get_wtime()

  self_it = self%it
  self_new = self%new
  self_time = self%time

  call openmp_work (end_time, self%nv, self%gn, self%mem, self_it, self_new, slope_type, csound, self_time, self%dtime, self%lb, self%ub, self%gamma, self%ds, self%g1, self%lo, self%uo, self%g2, self%li, self%ui)

  print *, "Finished solving...", omp_get_wtime() - time

  self%it    = self_it
  self%new   = self_new
  self%time  = self_time
  
#ifdef DBG 
  print '(f12.6,2x,a,i4,f8.2)', &
       omp_get_wtime()-start, 'solver_t%update: id, time =', self%id, self%time
#endif  
  !call self%counter_update
END SUBROUTINE update

!===============================================================================
!> Allocate work arrays
!===============================================================================
SUBROUTINE alloc (self_gn, self_nv, self_prim, self_left, self_rght, self_grad, self_flux)
  implicit none
  integer, dimension(3):: self_gn
  integer:: self_nv
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_left
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_rght
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_grad
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_flux
  !$omp declare target
  
  !class(solver_t):: self
  !integer, save:: itimer=0
  !.............................................................................
  !call tracebegin ('solver_t%alloc', itimer=itimer)
  !if (.not.associated(self%prim)) then
    allocate (self_prim (self_gn(1), self_gn(2), self_gn(3), self_nv   ))
    allocate (self_left (self_gn(1), self_gn(2) ,self_gn(3), self_nv   ))
    allocate (self_rght (self_gn(1), self_gn(2) ,self_gn(3), self_nv   ))
    allocate (self_grad (self_gn(1), self_gn(2) ,self_gn(3), self_nv, 3))
    allocate (self_flux (self_gn(1), self_gn(2) ,self_gn(3), self_nv, 3))
    !---------------------------------------------------------------------------
    ! Set to zero once, to allow undefined guard zone values to be used
    !---------------------------------------------------------------------------
    self_prim = 0.0
    self_left = 0.0
    self_rght = 0.0
    self_grad = 0.0
    self_flux = 0.0
  !end if
  !call traceend (itimer)
END SUBROUTINE alloc

!===============================================================================
!> Deallocate work arrays
!===============================================================================
SUBROUTINE dealloc (self_prim, self_left, self_rght, self_grad, self_flux)
  implicit none
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_left
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_rght
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_grad
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_flux
  !$omp declare target
  !class(solver_t):: self
  !integer, save:: itimer=0
  !.............................................................................
  deallocate (self_prim)
  deallocate (self_left)
  deallocate (self_rght)
  deallocate (self_grad)
  deallocate (self_flux)
END SUBROUTINE dealloc

!===============================================================================
!> Conservative to primitive transformation, using log(d), U, and Eint
!> Valid in full interval: [lb,ub] = [li-2,ui+2]
!===============================================================================
SUBROUTINE ctoprim (self_it, self_new, self_lb, self_ub, self_mem, self_prim, self_gamma)
  implicit none
  integer:: self_it, self_new
  integer, dimension(3):: self_lb, self_ub
  real, pointer:: self_mem(:,:,:,:,:,:)
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real(8):: self_gamma

  !$omp declare target
  !class(solver_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, l(3), u(3)
  !type(lock_t), pointer:: lock
  real:: u2, u2_max, c2, g2
  !integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  !call tracebegin ('solver_t%ctoprim', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Copy everything to new time slot, lock protected
  !-----------------------------------------------------------------------------

  ! TODO avoid copying self_mem in and out
  
  !lock => self%mem_lock(self%it)
  !call lock%set ('ctoprim')
  self_mem(:,:,:,:,self_new,1) = self_mem(:,:,:,:,self_it,1)
  !call lock%unset ('ctoprim')
  !-----------------------------------------------------------------------------
  l = self_lb
  u = self_ub

  !$omp teams distribute parallel do schedule(static,1) shared(self_prim, self_mem, l, u, self_new) private(iv, i3, i2) default(none) collapse(2)
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    self_prim(:,i2,i3,1) = log(self_mem(:,i2,i3, 1,self_new,1))
    do iv=2,5
      self_prim(:,i2,i3,iv)  = self_mem(:,i2,i3,iv,self_new,1) &
                        / self_mem(:,i2,i3, 1,self_new,1)
    end do
  end do
  end do
  !$omp end teams distribute parallel do

  !-----------------------------------------------------------------------------
  ! Subtract kinetic energy per unit mass and compute self%u_max
  !-----------------------------------------------------------------------------
  if (self_gamma == 1d0) then
    g2 = 1.0
  else
    g2 = (self_gamma-1d0)*self_gamma
  end if
  u2_max = 0.0
  
  !$omp teams distribute parallel do schedule(static,1) shared(self_prim, l, u) private(i1, i3, i2, u2) default(none) collapse(3)
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    do i1=l(1),u(1)
      u2 = self_prim(i1,i2,i3,2)**2 + &
           self_prim(i1,i2,i3,3)**2 + &
           self_prim(i1,i2,i3,4)**2
      self_prim(i1,i2,i3,5) = self_prim(i1,i2,i3,5) - 0.5*u2
      !c2 = self_prim(i1,i2,i3,5)*g2
      !u2_max = max (u2_max, u2+c2)
      self_prim(i1,i2,i3,5) = log(self_prim(i1,i2,i3,5))
    end do
  end do
  end do
  !$omp end teams distribute parallel do

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
SUBROUTINE slopes (self_nv, self_lb, self_ub, local_slope_type, self_prim, self_grad)
  implicit none
  !class(solver_t):: self
  integer:: self_nv
  integer, dimension(3):: self_lb, self_ub
  integer:: local_slope_type
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_grad
  !$omp declare target
  !.............................................................................
  integer:: i1, i2, i3, iv, dir, l(3), u(3)
  !integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  !call tracebegin ('solver_t%slopes', itimer=itimer)
  l = self_lb+1
  u = self_ub-1
  if (local_slope_type == 1) then
    call minmod(self_nv, self_prim, self_grad, l, u)
  else if (local_slope_type >= 2) then
    call moncen(self_nv, self_prim, self_grad, l, u)
  else
    call central(self_nv, self_grad, self_prim, l, u)
  end if
  !call traceend (itimer)
contains
!===============================================================================
subroutine central(self_nv, self_grad, self_prim, l, u)
  implicit none
  integer:: self_nv
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_grad
  integer:: l(3), u(3)
  !$omp declare target

  do iv=1,self_nv
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    do i1=l(1),u(1)
      self_grad (i1,i2,i3,iv,1) = 0.5*(self_prim(i1+1,i2,i3,iv) - self_prim(i1-1,i2,i3,iv))
      self_grad (i1,i2,i3,iv,2) = 0.5*(self_prim(i1,i2+1,i3,iv) - self_prim(i1,i2-1,i3,iv))
      self_grad (i1,i2,i3,iv,3) = 0.5*(self_prim(i1,i2,i3+1,iv) - self_prim(i1,i2,i3-1,iv))
    end do
  end do
  end do
  end do
end subroutine central
!===============================================================================
subroutine minmod(self_nv, self_prim, self_grad, l, u)
  implicit none
  integer:: self_nv
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_grad
  integer:: l(3), u(3)
  !$omp declare target

  real:: dlft, drgt ! private to each thread
  do iv=1,self_nv
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !$omp simd
    do i1=l(1),u(1)
      ! slopes in first coordinate direction
      dlft = self_prim(i1  ,i2,i3,iv) - self_prim(i1-1,i2,i3,iv)
      drgt = self_prim(i1+1,i2,i3,iv) - self_prim(i1  ,i2,i3,iv)
      self_grad(i1,i2,i3,iv,1) = merge(0.0,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=0.0)
      ! slopes in second coordinate direction
      dlft = self_prim(i1,i2  ,i3,iv) - self_prim(i1,i2-1,i3,iv)
      drgt = self_prim(i1,i2+1,i3,iv) - self_prim(i1,i2  ,i3,iv)
      self_grad(i1,i2,i3,iv,2) = merge(0.0,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=0.0)
      ! slopes in third coordinate direction
      dlft = self_prim(i1,i2,i3  ,iv) - self_prim(i1,i2,i3-1,iv)
      drgt = self_prim(i1,i2,i3+1,iv) - self_prim(i1,i2,i3  ,iv)
      self_grad(i1,i2,i3,iv,3) = merge(0.0,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=0.0)
    end do
  end do
  end do
  end do
end subroutine minmod
!===============================================================================
subroutine moncen(self_nv, self_prim, self_grad, l, u)
  implicit none
  integer:: self_nv
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_grad
  integer:: l(3), u(3)
  !$omp declare target

  real:: dlft,  drgt, dcen, dsgn, slop, dlim, c1, c2 ! private
  c2 = 0.5 ! shared
  !$omp teams distribute parallel do schedule(static,1) shared(self_prim, c2, self_grad, l, u, self_nv) private(iv, i3, i2, i1, dlft, drgt, dcen, dsgn, slop, dlim) default(none) collapse(4)
  do iv=1,self_nv
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    do i1 = l(1),u(1)
      ! slopes in first coordinate direction
      dlft = self_prim(i1  ,i2,i3,iv) - self_prim(i1-1,i2,i3,iv)
      drgt = self_prim(i1+1,i2,i3,iv) - self_prim(i1  ,i2,i3,iv)
      dcen = c2*(dlft+drgt)
      dsgn = sign(1.0, dcen)
      slop = min(abs(dlft),abs(drgt))
      dlim = merge(0.0,slop,(dlft*drgt)<=0.0)
      self_grad(i1,i2,i3,iv,1) = dsgn*min(dlim,abs(dcen))
      ! slopes in second coordinate direction
      dlft = self_prim(i1,i2  ,i3,iv) - self_prim(i1,i2-1,i3,iv)
      drgt = self_prim(i1,i2+1,i3,iv) - self_prim(i1,i2  ,i3,iv)
      dcen = c2*(dlft+drgt)
      dsgn = sign(1.0,dcen)
      slop = min(abs(dlft),abs(drgt))
      dlim = merge(0.0,slop,(dlft*drgt)<=0.0)
      self_grad(i1,i2,i3,iv,2) = dsgn*min(dlim,abs(dcen))
      ! slopes in third coordinate direction
      dlft = self_prim(i1,i2,i3  ,iv) - self_prim(i1,i2,i3-1,iv)
      drgt = self_prim(i1,i2,i3+1,iv) - self_prim(i1,i2,i3  ,iv)
      dcen = c2*(dlft+drgt)
      dsgn = sign(1.0,dcen)
      slop = min(abs(dlft),abs(drgt))
      dlim = merge(0.0,slop,(dlft*drgt)<=0.0)
      self_grad(i1,i2,i3,iv,3) = dsgn*min(dlim,abs(dcen))
    end do
  end do
  end do
 end do
 !$omp end teams distribute parallel do
end subroutine moncen
END SUBROUTINE slopes

!===============================================================================
!> Predictor step for primitive variables; should produce left and right states
!> in intervals [li,uo] for left and [lo,uo] for right, so that after a shift
!> of right one index, the right state is valid in the interval [lo,ui]
!===============================================================================
SUBROUTINE predict (self_g1, self_lo, self_uo, self_dtime, self_new, self_ds, self_prim, self_grad, self_mem)
  implicit none
  real(8):: self_g1
  integer:: self_lo(3), self_uo(3)
  real(8):: self_dtime
  integer:: self_new
  real(8):: self_ds(3)
  real, pointer:: self_mem(:,:,:,:,:,:)
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_grad  
  !$omp declare target
  !class(solver_t):: self
  real, pointer, contiguous, dimension(:,:,:):: &
    ln_1, ln_2, ln_3, &
    ee_1, ee_2, ee_3, &
    u1_1, u1_2, u1_3, &
    u2_1, u2_2, u2_3, &
    u3_1, u3_2, u3_3, &
    ln, ee, u1, u2, u3, d
  real:: divu, dtd1, dtd2, dtd3, g1, dt2
  real:: ff(size(self_prim,1),3)
  integer:: i1, i2, i3, l(3), u(3)
  !integer, save:: itimer=0
  !.............................................................................
  !call tracebegin ('solver_t%predict', itimer=itimer)
  ln_1 => self_grad(:,:,:,1,1)
  ln_2 => self_grad(:,:,:,1,2)
  ln_3 => self_grad(:,:,:,1,3)
  u1_1 => self_grad(:,:,:,2,1)
  u1_2 => self_grad(:,:,:,2,2)
  u1_3 => self_grad(:,:,:,2,3)
  u2_1 => self_grad(:,:,:,3,1)
  u2_2 => self_grad(:,:,:,3,2)
  u2_3 => self_grad(:,:,:,3,3)
  u3_1 => self_grad(:,:,:,4,1)
  u3_2 => self_grad(:,:,:,4,2)
  u3_3 => self_grad(:,:,:,4,3)
  ee_1 => self_grad(:,:,:,5,1)
  ee_2 => self_grad(:,:,:,5,2)
  ee_3 => self_grad(:,:,:,5,3)
  ln   => self_prim(:,:,:,1)
  u1   => self_prim(:,:,:,2)
  u2   => self_prim(:,:,:,3)
  u3   => self_prim(:,:,:,4)
  ee   => self_prim(:,:,:,5)
  d    => self_mem(:,:,:,1,self_new,1)
  !-----------------------------------------------------------------------------
  ! Evolution of primitive variables over half a time step, with
  ! \partial\ln(d) / \partial t = -u*grad(ln_d) - div(u)
  !     \partial u / \partial t = -u*grad(u) - grad(P)/d
  !                     grad(P) = d*grad(e) + e*grad(d)
  !                             = d*grad(e) + e*d*grad(\ln(d))
  !     \partial e / \partial t = -u*grad(e) - (P/d)*div(U)
  !-----------------------------------------------------------------------------
  dt2  = self_dtime*0.5
  dtd1 = dt2/self_ds(1)
  dtd2 = dt2/self_ds(2)
  dtd3 = dt2/self_ds(3)
  g1 = self_g1
  l = self_lo
  u = self_uo

  !$omp teams distribute parallel do schedule(static,1) shared(ff, l, u, g1, ee_2, ee_3, u3_2, u3_1, u2_3, u2_1, dt2, dtd1, dtd2, dtd3, ee, ee_1, u1_3, u1_2, ln_3, u3, ln_2, u2, ln_1, u1, ln, u3_3, u2_2, u1_1) private(i3, i2, i1, divu) default(none) collapse(2)
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    ff(:,:) = 0.0
    !if (allocated(self_force_per_unit_mass)) then
    !  ff(:,1) = ff(:,1) + self_force_per_unit_mass(:,i2,i3,1)
    !  ff(:,2) = ff(:,2) + self_force_per_unit_mass(:,i2,i3,2)
    !  ff(:,3) = ff(:,3) + self_force_per_unit_mass(:,i2,i3,3)
    !end if
    !if (allocated(self_force_per_unit_volume)) then
    !  ff(:,1) = ff(:,1) + self_force_per_unit_volume(:,i2,i3,1)/d(:,i2,i3)
    !  ff(:,2) = ff(:,2) + self_force_per_unit_volume(:,i2,i3,2)/d(:,i2,i3)
    !  ff(:,3) = ff(:,3) + self_force_per_unit_volume(:,i2,i3,3)/d(:,i2,i3)
    !end if
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
  !$omp end teams distribute parallel do  

  !call traceend (itimer)
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
SUBROUTINE riemn3d (self_g1, self_g2, self_lo, self_uo, self_nv, self_left, self_rght, self_prim, self_grad, self_flux)
  implicit none
  !class(solver_t):: self
  real(8):: self_g1, self_g2
  integer, dimension(3):: self_lo, self_uo
  integer:: self_nv
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_left
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_rght
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_grad
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_flux
  !$omp declare target
  !.............................................................................
  integer:: l(3), u(3), i1, i2, i3, iv
!  integer, save:: itimer=0, itimer1=0
  !-----------------------------------------------------------------------------
!  call tracebegin ('solver_t%riemn3d', itimer=itimer)
  l = self_lo
  u = self_uo
  !-----------------------------------------------------------------------------
  ! X-direction
  !-----------------------------------------------------------------------------

  !$omp teams distribute parallel do schedule(static,1) shared(self_grad, self_prim, self_left, self_rght, self_nv, l, u) private(iv, i3, i2, i1) default(none) collapse(4)
  do iv=1,self_nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    self_left(i1,i2,i3,iv) = self_prim(i1,i2,i3,iv) + self_grad(i1,i2,i3,iv,1)*0.5
    self_rght(i1,i2,i3,iv) = self_prim(i1,i2,i3,iv) - self_grad(i1,i2,i3,iv,1)*0.5
  end do; end do; end do; end do
  !$omp end teams distribute parallel do 
!                    call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  call solve (self_g1, self_g2, 1, l+1, u, self_left, self_rght, self_flux(:,:,:,:,1))
!                    call traceend (itimer1)
  !-----------------------------------------------------------------------------
  ! Y-direction
  !-----------------------------------------------------------------------------

  !$omp teams distribute parallel do schedule(static,1) shared(self_grad, self_prim, self_left, self_rght, self_nv, l, u) private(iv, i3, i2, i1) default(none) collapse(4)
  do iv=1,self_nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    self_left(i1,i2,i3,iv) = self_prim(i1,i2,i3,iv) + self_grad(i1,i2,i3,iv,2)*0.5
    self_rght(i1,i2,i3,iv) = self_prim(i1,i2,i3,iv) - self_grad(i1,i2,i3,iv,2)*0.5
  end do; end do; end do; end do
  !$omp end teams distribute parallel do 
!                    call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  call solve (self_g1, self_g2, 2, l+1, u, self_left, self_rght, self_flux(:,:,:,:,2))
!                    call traceend (itimer1)
  !-----------------------------------------------------------------------------
  ! Z-direction
  !-----------------------------------------------------------------------------

  !$omp teams distribute parallel do schedule(static,1) shared(self_grad, self_prim, self_left, self_rght, self_nv, l, u) private(iv, i3, i2, i1) default(none) collapse(4)
  do iv=1,self_nv; do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
    self_left(i1,i2,i3,iv) = self_prim(i1,i2,i3,iv) + self_grad(i1,i2,i3,iv,3)*0.5
    self_rght(i1,i2,i3,iv) = self_prim(i1,i2,i3,iv) - self_grad(i1,i2,i3,iv,3)*0.5
  end do; end do; end do; end do
  !$omp end teams distribute parallel do 
!                    call tracebegin ('riemann_t%'//trim(method), itimer=itimer1)
  call solve (self_g1, self_g2, 3, l+1, u, self_left, self_rght, self_flux(:,:,:,:,3))
!                    call traceend (itimer1)
  !-----------------------------------------------------------------------------
!  call traceend (itimer)
END SUBROUTINE riemn3d

!===============================================================================
!> Update conservative variables with flux divergences, where flux(i1,:,:,:,1)
!> is centered at i1-0.5, so needs to be defined in the interval [li,uo]
!===============================================================================
SUBROUTINE divflux (self_new, self_dtime, self_ds, self_li, self_ui, self_nv, self_mem, self_flux)
  implicit none
  real(8):: self_dtime
  real(8):: self_ds(3)
  integer, dimension(3):: self_li, self_ui
  integer:: self_nv, self_new
  real, pointer:: self_mem(:,:,:,:,:,:)
  real, dimension(:,:,:,:,:), pointer, contiguous:: self_flux
  !$omp declare target
  !class(solver_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, dir, l(3), u(3)
  real:: dtd1, dtd2, dtd3
  !integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  !call tracebegin ('solver_t%divflux', itimer=itimer)
  dtd1 = self_dtime/self_ds(1)
  dtd2 = self_dtime/self_ds(2)
  dtd3 = self_dtime/self_ds(3)
  l = self_li
  u = self_ui

  !$omp teams distribute parallel do schedule(static,1) shared(self_mem, self_new, self_nv, l, u, self_flux, dtd1, dtd2, dtd3) private(iv, i3, i2, i1) default(none) collapse(4)
  do iv=1,self_nv
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    self_mem(i1,i2,i3,iv,self_new,1) = &
    self_mem(i1,i2,i3,iv,self_new,1)   &
    - dtd1*(self_flux(i1+1,i2  ,i3  ,iv,1) - self_flux(i1,i2,i3,iv,1)) &
    - dtd2*(self_flux(i1  ,i2+1,i3  ,iv,2) - self_flux(i1,i2,i3,iv,2)) &
    - dtd3*(self_flux(i1  ,i2  ,i3+1,iv,3) - self_flux(i1,i2,i3,iv,3))
  end do
  end do
  end do
  end do
  !$omp end teams distribute parallel do
  !call traceend (itimer)
END SUBROUTINE divflux

!===============================================================================
!> Add source terms
!===============================================================================
SUBROUTINE sources (local_csound, self_new, self_gamma, self_mem, self_prim)
  implicit none
  real:: local_csound
  integer:: self_new
  real(8):: self_gamma
  real, pointer:: self_mem(:,:,:,:,:,:)
  real, dimension(:,:,:,:)  , pointer, contiguous:: self_prim
  !$omp declare target
  !class(solver_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iv, dir
  real:: dtd1, dtd2, dtd3
  real, dimension(:,:,:), pointer, contiguous:: p1, p2, p3, d, e
  logical:: forcing
  !integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  !call tracebegin ('solver_t%sources', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! If adjusting total energy, save kinetic part in self%prim(:,:,:,5)
  !-----------------------------------------------------------------------------
  d  => self_mem(:,:,:,1,self_new,1)
  p1 => self_mem(:,:,:,2,self_new,1)
  p2 => self_mem(:,:,:,3,self_new,1)
  p3 => self_mem(:,:,:,4,self_new,1)
  e  => self_mem(:,:,:,5,self_new,1)
  forcing = .false.
  !forcing = allocated(self_force_per_unit_mass) .or. &
  !!          allocated(self%force_per_unit_volume)
  if (self_gamma /= 1d0 .and. forcing) then
    self_prim(:,:,:,5) = p1**2 + p2**2 + p3**2
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
  if (self_gamma == 1d0) then
    e = d*local_csound**2 + 0.5*(p1**2 + p2**2 + p3**2)/d
  else if (forcing) then
    e = e + 0.5*(p1**2 + p2**2 + p3**2 - self_prim(:,:,:,5))/d
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
