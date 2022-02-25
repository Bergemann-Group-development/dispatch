!===============================================================================
!> HLLD solver for GPUs (and CPUs).  To avoid potential issues with derived
!> data types, modules are (for now) only used to hold data and procedures
!>
!> call sequence inside update():
!>   ctoprim     :  compute primitive variables
!>   slopes      :  compute slopes of primitive variables
!>   preict      :  compute time derivatives and half-step primitives
!>   solve_1d    :  compute face values and solve for HD fluxes
!>   div_flux    :  update conserved variables with fluxes
!>
!> To minimize memory motions (also internally on the GPU), the left and right
!> face value arrays are 2-dimensional, with running index over the 1st dim,
!> which allows the 1-D solver to parallelize over the two redundant dimensions,
!> based only on access to the primitive variables and their (half) slopes.
!>
!> The results is walking through memory 7 times:
!> 1. ctoprim     uold -> prim
!> 2. slopes      prim -> dprm
!> 3. predict     prim -> prim
!> 4. solve_1d    prim,dprim -> flux
!> 6. divflux     flux,uold  -> unew
!===============================================================================
MODULE mhd_mod
  USE io_mod
  USE os_mod
  USE trace_mod
  USE omp_lock_mod
  USE extras_mod
  USE eos_mod
  USE index_mod
  USE scratch_mod
  ! USE divb_cleaner_mod
  implicit none
  private
  type, public, extends(extras_t):: mhd_t
    real, pointer, contiguous, dimension(:,:,:,:)   :: prim
    real, pointer, contiguous, dimension(:,:,:,:,:) :: dprim, flux
    real, pointer, contiguous, dimension(:,:,:)     :: face
    real                                            :: cdtd = 0.5
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
    procedure:: post_update
    procedure:: update
    procedure, nopass:: ctoprim
    procedure:: div_flux
  end type
  integer:: r1=1, u1=2, u2=3, u3=4, e1=5
  real, parameter:: one=1d0, half=0.5, zero=0.0
#ifdef OFFLOAD
  !$omp declare target(n_ctoprim,n_solve_1d,n_solve_2d,n_predict,n_slopes)
  !$omp declare target(slope_type)
#else
  !save:: prim, gpu1, gpu2, gpu3, gpu4, gpu5
  !!omp threadprivate (prim,gpu1,gpu2,gpu3,gpu5)
#endif
  integer :: verbose        = 0
  integer :: unit           = 0
  integer :: courant_type   = 0
  real    :: courant_d      = 1.0
  real    :: slope_type     = 3.5 
  real    :: smallr         = 0.0
  real    :: smallc         = 0.1
  real    :: courant        = 0.25
  real    :: csound         = 1.00
  real(8) :: gamma          = 1.0
  logical :: unsigned       = .true.
  logical :: detailed_timer = .false.
  logical, save:: isothermal =.false.
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE pre_init (self)
  class(mhd_t), target:: self
  integer:: iostat
  namelist /an_mhd_params/ verbose, courant, courant_type, slope_type, csound, &
                           courant_d, smallr, smallc, unsigned, detailed_timer, isothermal
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, an_mhd_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, an_mhd_params)
    end if
    gamma = eos%gamma
  end if
  call io%lock%unset
  self%courant = courant
  self%ng = 2
  self%nt = 5
  self%nv = 5
  self%nw = 1
  self%kind  = 'AN_hllc'
  self%type  = 'mhd_t'
  self%etype = 'thermal'
  self%staggered = .false.
  call self%idx%init (energy=5, mhd=.false.)
  e1 = self%idx%e
  call self%extras_t%pre_init ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize
!===============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  integer:: i
  !-----------------------------------------------------------------------------
  call self%extras_t%init
  do i=1,3
    self%mesh(i)%h = 0d0
  end do
  self%unsigned(self%idx%d) = unsigned
  self%unsigned(self%idx%e) = unsigned
  self%pervolume(self%idx%px) = unsigned
  self%pervolume(self%idx%py) = unsigned
  self%pervolume(self%idx%pz) = unsigned
END SUBROUTINE init

!===============================================================================
!> Allocate scratch
!===============================================================================
SUBROUTINE pre_update (self)
  class(mhd_t), target:: self
  integer:: m(6)
  integer:: l(4), u(4)
  real, contiguous, pointer, dimension(:,:,:,:):: uold, unew
  real:: u_max
  !-----------------------------------------------------------------------------
  call self%extras_t%pre_update ()
  m = shape(self%mem)
  m(4) = self%nv
  call os%allocate (self%dprim, [m(1),m(4),3,m(2),m(3)])
  call os%allocate (self%flux , [m(1),m(2),m(3),m(4),3])
  call os%allocate (self%prim , [m(1),m(4),  m(2),m(3)])
  call os%allocate (self%face , [m(1),m(4),2])
  !-----------------------------------------------------------------------------
  uold => self%mem(:,:,:,:,self%it,1)
  unew => self%mem(:,:,:,:,self%new,1)
  l(4) = 1
  u(4) = self%nv
  !-----------------------------------------------------------------------------
  ! The primitive values are needed in one cell more than the slopes, and to
  ! compute them requires access to one more cell below, so from li-3
  !-----------------------------------------------------------------------------
  l(1:3) = self%mesh%lb
  u(1:3) = self%mesh%ub
  call ctoprim  (self, l, u, uold, self%prim, gamma, u_max)
  self%u_max = u_max
  !-----------------------------------------------------------------------------
  ! Pass the call on to extras
  !-----------------------------------------------------------------------------
  call self%extras_t%courant_condition ()
  call self%gpatch_t%pre_update ()
  call self%extras_t%pre_update ()
END SUBROUTINE pre_update

!===============================================================================
!> Update organizer / public entry point
!>
!> In the sketch of guard zones below, 1-2 are the guard cells, + are internal
!> 1---2---+................+---2---1        ! full range
!> 1---2---+................+---2---1        ! ctoprim valid
!>     2---+................+---2            ! slopes valid
!>     2---+................+---2            ! right valid
!>         +................+---2            ! left valid
!>         +................+---2            ! riemann valid
!>         +................+                ! div(flux) valid
!===============================================================================
SUBROUTINE update (self)
#ifdef OFFLOAD
  !$omp declare target
#endif
  class(mhd_t)              :: self
  integer                   :: l(4), u(4), j
  real, contiguous, pointer :: uold(:,:,:,:), unew(:,:,:,:)
  type(lock_t),     pointer :: lock
  real                      :: dtds, dt2, dt
  integer,             save :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%update', itimer=itimer)
  call self%output ()
  l(1:3) = self%mesh%lb+1
  u(1:3) = self%mesh%ub-1
  l(4) = 1
  u(4) = self%nv
  call slopes   (l, u, self%prim, self%dprim)
  !-----------------------------------------------------------------------------
  ! The predictor needs to compute values one cell lower than the flux
  !-----------------------------------------------------------------------------
  l(1:3) = self%mesh%lb+1
  u(1:3) = self%mesh%ub-1
  dtds = self%dtime/self%ds(1)
  dt2 = self%dtime/2.0
  call predict  (l, u, self%prim, self%dprim, gamma, dtds, dt2)
  !-----------------------------------------------------------------------------
  ! The fluxes need to be computed in the interval li-uo, to allow dd_up
  !-----------------------------------------------------------------------------
  l(1:3) = self%mesh%lb+2
  u(1:3) = self%mesh%ub-1
  call solve_1d (l, u, self%prim, self%dprim, self%flux, self%face, gamma)
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Post update, allowing extras_t to intervene
!===============================================================================
SUBROUTINE post_update (self)
  class(mhd_t)              :: self
  integer                   :: l(4), u(4)
  real, contiguous, pointer :: uold(:,:,:,:), unew(:,:,:,:)
  type(lock_t),     pointer :: lock
  real                      :: dtds, dt
  integer,             save :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%post_update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Allow extras%post_update to modify fluxes and EMFs
  !-----------------------------------------------------------------------------
  call self%extras_t%post_update ()
  !-----------------------------------------------------------------------------
  ! We need the updated variables in the interval li-ui
  !-----------------------------------------------------------------------------
  l(1:3) = self%mesh%li
  u(1:3) = self%mesh%ui
  l(4) = 1
  u(4) = self%nv
  !-----------------------------------------------------------------------------
  ! Divergence of fluxes and curl of EMFs
  !-----------------------------------------------------------------------------
  uold => self%mem(:,:,:,:,self%it,1)
  unew => self%mem(:,:,:,:,self%new,1)
  dtds = self%dtime/self%ds(1)
  lock => self%mem_lock(self%new)
                                                call lock%set ('mhd_t%div_flux')
  call div_flux (self, l, u, uold, unew, self%flux, dtds)
  !-----------------------------------------------------------------------------
  ! External sources
  !-----------------------------------------------------------------------------
  dt = self%dtime
  call sources  (l, u, unew, dt)
  !-----------------------------------------------------------------------------
  ! Handle isothermal cases
  !-----------------------------------------------------------------------------
  if (eos%gamma == 1d0) then
    unew(l(1):u(1),l(2):u(2),l(3):u(3),e1) = &
    unew(l(1):u(1),l(2):u(2),l(3):u(3),r1) * csound**2
  end if
                                              call lock%unset ('mhd_t%div_flux')
  !-----------------------------------------------------------------------------
  ! Wind up
  !-----------------------------------------------------------------------------
  call os%deallocate (self%prim)
  call os%deallocate (self%dprim)
  call os%deallocate (self%flux)
  call os%deallocate (self%face)
  call self%extras_t%post_update ()
  call self%gpatch_t%counter_update ()
  !.............................................................................
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE post_update

!===============================================================================
!> Conservative to primitive variables. Like all procedures below, it does not
!> refer to any structure components -- information is passed via argument lists
!===============================================================================
SUBROUTINE ctoprim (self, l, u, uold, prim, gamma, u_max)
#ifdef OFFLOAD
  !$omp declare target
#endif
  class(mhd_t)     :: self
  integer          :: l(4), u(4)
  real, contiguous :: uold(:,:,:,:), prim(:,:,:,:)
  real(8)          :: gamma
  real             :: u_max, u_fast, u_gas_1, u_gas_2, cs2 
  !.............................................................................
  integer:: i1, i2, i3, n_smalld, n_smallp
  real:: d, e, px, py, pz, ux, uy, uz, g0, g1, pg, u2_max, pg0, pg1, pg0_min
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%ctoprim', itimer=itimer)
  g0 = gamma
  g1 = merge (1d0, gamma-1d0, gamma==1d0)
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1) &
  !$omp private(i1,i2,i3,d,e,px,py,pz)
#endif
  u_max = 0.0
  u2_max = 0.0
  n_smalld = 0
  n_smallp = 0
  cs2 = csound**2
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !$omp simd
    do i1=l(1),u(1)
      d  = uold(i1,i2,i3,r1)
      e  = uold(i1,i2,i3,e1)
      px = uold(i1,i2,i3,u1)
      py = uold(i1,i2,i3,u2)
      pz = uold(i1,i2,i3,u3)
      n_smalld = merge (n_smalld+1, n_smalld, d<smallr)
      d = max(d,smallr)
      prim(i1,r1,i2,i3) = log(d)
      ux = px/d
      uy = py/d
      uz = pz/d
      prim(i1,u1,i2,i3) = ux
      prim(i1,u2,i2,i3) = uy
      prim(i1,u3,i2,i3) = uz
      pg0 = g1*e
      pg0_min = smallc**2/g0*d**g0
      n_smallp = merge (n_smallp+1, n_smallp, pg0<pg0_min)
      pg0 = max(pg0,pg0_min)
      pg1 = cs2*d
      pg = merge (pg1, pg0, gamma==1d0)
      u_fast = sqrt((gamma*pg)/d)
      u_gas_1 = max(abs(ux),abs(uy),abs(uz))
      u_gas_2 = (abs(ux)+abs(uy)+abs(uz))*0.3333333
      u_max = max(u_fast + merge (u_gas_1, u_gas_2, courant_type==1), u_max)
      prim(i1,e1,i2,i3) = log(pg)
    end do
  end do
  end do
  !.............................................................................
  if (verbose > 0) then
    unit = merge (stdout, stdlog, verbose > 1)
    if (n_smalld > 0) &
      write (unit,'(a,2i6,1p,g14.5)') &
        'mhd_t%div_flux: id,n_smalld,time=', self%id, n_smalld, self%time
    if (n_smallp > 0) &
      write (unit,'(a,2i6,1p,g14.5)') &
        'mhd_t%div_flux: id,n_smallp,time=', self%id, n_smallp, self%time
  end if
  !.............................................................................
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE ctoprim

!===============================================================================
!> 3-D positivity preserving slope limiter
!===============================================================================
SUBROUTINE slopes (l, u, q, dq)
#ifdef OFFLOAD
  !$omp declare target
#endif
  integer:: l(4), u(4)
  real:: q(:,:,:,:), dq(:,:,:,:,:)
  real::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real::vmin,vmax,dfx,dfy,dfz,dff,idff,dlim
  integer:: i, j, k, n
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%slopes', itimer=itimer)
#ifdef OFFLOAD
  !$omp parallel do collapse(3) schedule(static,1)
#endif
  do n=l(4),u(4)
  do k=l(3),u(3)
  do j=l(2),u(2)
    !$omp simd
    do i=l(1),u(1)
      vmin = min(q(i-1,n,j-1,k-1), q(i,n,j-1,k-1), q(i+1,n,j-1,k-1), &
                 q(i-1,n,j  ,k-1), q(i,n,j  ,k-1), q(i+1,n,j  ,k-1), &
                 q(i-1,n,j+1,k-1), q(i,n,j+1,k-1), q(i+1,n,j+1,k-1), &
                 q(i-1,n,j-1,k  ), q(i,n,j-1,k  ), q(i+1,n,j-1,k  ), &
                 q(i-1,n,j  ,k  ), q(i,n,j  ,k  ), q(i+1,n,j  ,k  ), &
                 q(i-1,n,j+1,k  ), q(i,n,j+1,k  ), q(i+1,n,j+1,k  ), &
                 q(i-1,n,j-1,k+1), q(i,n,j-1,k+1), q(i+1,n,j-1,k+1), &
                 q(i-1,n,j  ,k+1), q(i,n,j  ,k+1), q(i+1,n,j  ,k+1), &
                 q(i-1,n,j+1,k+1), q(i,n,j+1,k+1), q(i+1,n,j+1,k+1)) &
                -q(i,n,j,k)
      vmax = max(q(i-1,n,j-1,k-1), q(i,n,j-1,k-1), q(i+1,n,j-1,k-1), &
                 q(i-1,n,j  ,k-1), q(i,n,j  ,k-1), q(i+1,n,j  ,k-1), &
                 q(i-1,n,j+1,k-1), q(i,n,j+1,k-1), q(i+1,n,j+1,k-1), &
                 q(i-1,n,j-1,k  ), q(i,n,j-1,k  ), q(i+1,n,j-1,k  ), &
                 q(i-1,n,j  ,k  ), q(i,n,j  ,k  ), q(i+1,n,j  ,k  ), &
                 q(i-1,n,j+1,k  ), q(i,n,j+1,k  ), q(i+1,n,j+1,k  ), &
                 q(i-1,n,j-1,k+1), q(i,n,j-1,k+1), q(i+1,n,j-1,k+1), &
                 q(i-1,n,j  ,k+1), q(i,n,j  ,k+1), q(i+1,n,j  ,k+1), &
                 q(i-1,n,j+1,k+1), q(i,n,j+1,k+1), q(i+1,n,j+1,k+1)) &
                -q(i,n,j,k)
      dfx  = half*(q(i+1,n,j,k)-q(i-1,n,j,k))
      dfy  = half*(q(i,n,j+1,k)-q(i,n,j-1,k))
      dfz  = half*(q(i,n,j,k+1)-q(i,n,j,k-1))
      idff = (slope_type-2.0)/sqrt(dfx**2+dfy**2+dfz**2+tiny(1.0))
      !-------------------------------------------------------------------------
      ! For constant eg x-slope the expression below becomes 0.5*dfx unless
      ! one or more of the 26 nbor point differ by less than 1./(slope_type-2)
      ! from that.  So, for the std value 3.5 and slopes mainly along one
      ! direction, the reduction happens at 2/3 of the general slope
      !-------------------------------------------------------------------------
      dlim = 0.5*merge (min (1.0, min(abs(vmin),abs(vmax))*idff), 1.0, &
                        idff < 0.001*huge(1.0))
      dq(i,n,1,j,k) = dlim*dfx
      dq(i,n,2,j,k) = dlim*dfy
      dq(i,n,3,j,k) = dlim*dfz
    end do
  end do
  end do
  end do
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE slopes

!===============================================================================
!> Update the primitive variables half a timestep
!===============================================================================
SUBROUTINE predict (l, u, q, dq, gamma, dtds, dt2)
#ifdef OFFLOAD
  !$omp declare target
#endif
  integer                    :: l(4), u(4)
  real, dimension(:,:,:,:)   :: q      ! cell center
  real, dimension(:,:,:,:,:) :: dq     ! half slopes
  real(8)                    :: gamma  ! gas gamma
  real                       :: dtds   ! dt/ds
  real                       :: dt2    ! dt/2
  intent(in)                 :: dq
  contiguous                 :: q, dq
  !.............................................................................
  integer :: i1, i2, i3, i4
  real    :: g0, c
  real    :: dld, dlp, dux, duy, duz
  real    :: ld, ld_x, ld_y, ld_z, div_u
  real    :: lp, lp_x, lp_y, lp_z, tt
  real    :: ux, ux_x, ux_y, ux_z
  real    :: uy, uy_x, uy_y, uy_z
  real    :: uz, uz_x, uz_y, uz_z
  real    :: bx,       bx_y, bx_z, ex
  real    :: by, by_x,       by_z, ey
  real    :: bz, bz_x, bz_y      , ez
  real    :: ff(l(1):u(1),3), qq(l(1):u(1))
  integer :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%predict', itimer=itimer)
  g0 = gamma
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !---------------------------------------------------------------------------
    ! Include external forces in the predictor
    !---------------------------------------------------------------------------
    ff(:,:) = 0.0
    if (associated(scratch%force_per_unit_mass)) then
      !dir$ simd
      do i1=l(1),u(1)
        ff(i1,1) = ff(i1,1) + scratch%force_per_unit_mass(i1,i2,i3,1)*dt2
        ff(i1,2) = ff(i1,2) + scratch%force_per_unit_mass(i1,i2,i3,2)*dt2
        ff(i1,3) = ff(i1,3) + scratch%force_per_unit_mass(i1,i2,i3,3)*dt2
      end do
    end if
    if (associated(scratch%force_per_unit_volume)) then
      !dir$ simd
      do i1=l(1),u(1)
        c = exp(-q(i1,r1,i2,i3))
        ff(i1,1) = ff(i1,1) + scratch%force_per_unit_volume(i1,i2,i3,1)*dt2*c
        ff(i1,2) = ff(i1,2) + scratch%force_per_unit_volume(i1,i2,i3,2)*dt2*c
        ff(i1,3) = ff(i1,3) + scratch%force_per_unit_volume(i1,i2,i3,3)*dt2*c
      end do
    end if
    qq(:) = 0.0
    if (associated(scratch%heating_per_unit_mass)) then
      !dir$ simd
      do i1=l(1),u(1)
        qq(i1) = qq(i1) + scratch%heating_per_unit_mass(i1,i2,i3)*dt2
      end do
    end if
    if (associated(scratch%heating_per_unit_volume)) then
      !dir$ simd
      do i1=l(1),u(1)
        c = exp(-q(i1,r1,i2,i3))
        qq(i1) = qq(i1) + scratch%heating_per_unit_volume(i1,i2,i3)*dt2*c
      end do
    end if
    !$omp simd
    do i1=l(1),u(1)
      ld   =  q(i1,r1,i2,i3)
      lp   =  q(i1,e1,i2,i3)
      ux   =  q(i1,u1,i2,i3)
      uy   =  q(i1,u2,i2,i3)
      uz   =  q(i1,u3,i2,i3)
      !-------------------------------------------------------------------------
      ! Slope values are assumed to be changes across half a cell
      !-------------------------------------------------------------------------
      ld_x = dq(i1,r1,1,i2,i3)
      ld_y = dq(i1,r1,2,i2,i3)
      ld_z = dq(i1,r1,3,i2,i3)
      lp_x = dq(i1,e1,1,i2,i3)
      lp_y = dq(i1,e1,2,i2,i3)
      lp_z = dq(i1,e1,3,i2,i3)
      ux_x = dq(i1,u1,1,i2,i3)
      ux_y = dq(i1,u1,2,i2,i3)
      ux_z = dq(i1,u1,3,i2,i3)
      uy_x = dq(i1,u2,1,i2,i3)
      uy_y = dq(i1,u2,2,i2,i3)
      uy_z = dq(i1,u2,3,i2,i3)
      uz_x = dq(i1,u3,1,i2,i3)
      uz_y = dq(i1,u3,2,i2,i3)
      uz_z = dq(i1,u3,3,i2,i3)
      !-------------------------------------------------------------------------
      ! log density and log pressure evolution
      !-------------------------------------------------------------------------
      div_u = ux_x + uy_y + uz_z
      dld = - dtds*(ux*ld_x + uy*ld_y + uz*ld_z +    div_u)
      dlp = - dtds*(ux*lp_x + uy*lp_y + uz*lp_z + g0*div_u) 
      !-------------------------------------------------------------------------
      ! Velocity evolution
      !-------------------------------------------------------------------------
      tt = exp(lp-ld)
      dux = - dtds*(ux*ux_x + uy*ux_y + uz*ux_z + tt*lp_x) + ff(i1,1)
      duy = - dtds*(ux*uy_x + uy*uy_y + uz*uy_z + tt*lp_y) + ff(i1,2)
      duz = - dtds*(ux*uz_x + uy*uz_y + uz*uz_z + tt*lp_z) + ff(i1,3)
      !-------------------------------------------------------------------------
      ! Store results
      !-------------------------------------------------------------------------
      q(i1,r1,i2,i3) = ld + dld
      q(i1,e1,i2,i3) = lp + dlp
      q(i1,u1,i2,i3) = ux + dux
      q(i1,u2,i2,i3) = uy + duy
      q(i1,u3,i2,i3) = uz + duz
    end do
  end do
  end do
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE predict

!===============================================================================
!> Produce the face value that may be used to loop over calls to the Riemann-1D
!> solver.  Input arrays:
!>
!>   q: primitive variables
!>  dq: half slopes
!===============================================================================
SUBROUTINE solve_1d (l, u, q, dq, flux, a, gamma)
#ifdef OFFLOAD
  !$omp declare target
#endif
  integer          :: l(4), u(4)
  real, contiguous :: q(:,:,:,:), dq(:,:,:,:,:), flux(:,:,:,:,:), a(:,:,:)
  intent(in)       :: q, dq
  real(8)          :: gamma
  !.............................................................................
  integer :: i1, i2, i3, i4, j, j1, j2, j3, jj(3), iv
  integer :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%solve_1d', itimer=itimer)
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    do j=1,3
      jj(:) = 0
      jj(j) = 1
      j1 = jj(1)
      j2 = jj(2)
      j3 = jj(3)
      do i4=1,5
        !$omp simd
        do i1=l(1),u(1)
          a(i1,i4,1) =  q(i1-j1,i4,  i2-j2,i3-j3)  & ! qleft
                     + dq(i1-j1,i4,j,i2-j2,i3-j3)    ! ds/2 increment
          a(i1,i4,2) =  q(i1   ,i4,  i2   ,i3   )  & ! qrght
                     - dq(i1   ,i4,j,i2   ,i3   )    ! ds/2 increment
        end do
      end do
      !$omp simd
      do i1=l(1),u(1)
        a(i1,r1,1) = exp(a(i1,r1,1))
        a(i1,r1,2) = exp(a(i1,r1,2))
        a(i1,e1,1) = exp(a(i1,e1,1))
        a(i1,e1,2) = exp(a(i1,e1,2))
      end do
      call hllc (l, u, a(:,:,1), a(:,:,2), flux(:,i2,i3,:,:), j)
    end do
    end do
  end do
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE solve_1d

!===============================================================================
!> The RAMSES 1-D HLLD solver, verbatim.
!===============================================================================
SUBROUTINE hllc (lb, ub, qleft, qrght, flux, j)
#ifdef OFFLOAD
  !$omp declare target
#endif
  integer                               :: lb(4), ub(4)
  real, dimension(:,:),   intent(inout) :: qleft
  real, dimension(:,:),   intent(in   ) :: qrght
  real, dimension(:,:,:), intent(out  ) :: flux
  integer                               :: j
  real(8)                               :: gamma
  !.............................................................................
  real::sl,sr
  real::entho
  real::rl,pl,ul,vl,wl,ecinl,etotl,ptotl
  real::rr,pr,ur,vr,wr,ecinr,etotr,ptotr
  real::cfastl,rcl,rstarl
  real::cfastr,rcr,rstarr
  real::etotstarl,etotstarr
  real::ustar,ptotstar
  real::ro,uo,ptoto,etoto
  real::smallp
  integer ::ivar, i, icase, j1, j2, j3, v1, v2, v3
  logical :: error,detailed_timer
  !
  j1 = j
  j2 = 1 + mod(j1,3)
  j3 = 1 + mod(j2,3)
  v1 = u1+j1-1
  v2 = u1+j2-1
  v3 = u1+j3-1
  ! constants
  smallp = smallc**2/eos%gamma
  if (isothermal .or. eos%gamma==1d0) then
    entho=one
  else
    entho=one/(eos%gamma-one)
  end if
  error = .false.
  !
  do i=lb(1),ub(1)
    !
    ! Left & right variables
    rl=max(qleft(i,r1),smallr)
    rr=max(qrght(i,r1),smallr)
    Pl=max(qleft(i,e1),rl*smallp)
    Pr=max(qrght(i,e1),rr*smallp)
    ul=    qleft(i,v1)
    ur=    qrght(i,v1)
    vl=    qleft(i,v2)
    vr=    qrght(i,v2)
    wl=    qleft(i,v3)
    wr=    qrght(i,v3)
    !
    ecinl = half*rl*(ul*ul+vl*vl+wl*wl)
    ecinr = half*rr*(ur*ur+vr*vr+wr*wr)
    etotl = Pl*entho+ecinl
    etotr = Pr*entho+ecinr
    Ptotl = Pl
    Ptotr = Pr
    !
    ! Find the largest eigenvalues in the normal direction to the interface
    cfastl=sqrt(max(eos%gamma*Pl/rl,smallc**2))
    cfastr=sqrt(max(eos%gamma*Pr/rr,smallc**2))
    !
    ! Compute HLL wave speed
    SL=min(ul,ur)-max(cfastl,cfastr)
    SR=max(ul,ur)+max(cfastl,cfastr)
    !
    ! Compute lagrangian sound speed
    rcl=rl*(ul-SL)
    rcr=rr*(SR-ur)
    !
    ! Compute acoustic star state
    ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
    Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)
    !
    ! Left star region variables
    rstarl=rl*(SL-ul)/(SL-ustar)
    etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar)/(SL-ustar)
    !
    ! Right star region variables
    rstarr=rr*(SR-ur)/(SR-ustar)
    etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar)/(SR-ustar)
    !
    ! Sample the solution at x/t=0
!define MERGE
#ifdef MERGE
    ! -- default values
    ro=rr
    uo=ur
    Ptoto=Ptotr
    etoto=etotr
    ! -- 3rd if block
    ro   =merge(rstarr   ,ro     ,SR>0d0)
    uo   =merge(ustar    ,uo     ,SR>0d0)
    Ptoto=merge(Ptotstar ,Ptoto  ,SR>0d0)
    etoto=merge(etotstarr,etoto  ,SR>0d0)
    ! -- 2nd if block
    ro   =merge(rstarl   ,ro     ,ustar>0d0)
    uo   =merge(ustar    ,uo     ,ustar>0d0)
    Ptoto=merge(Ptotstar ,Ptoto  ,ustar>0d0)
    etoto=merge(etotstarl,etoto  ,ustar>0d0)
    ! -- 1sr if block
    ro   =merge(rl       ,ro     ,SL>0d0)
    uo   =merge(ul       ,uo     ,SL>0d0)
    Ptoto=merge(Ptotl    ,Ptoto  ,SL>0d0)
    etoto=merge(etotl    ,etoto  ,SL>0d0)
#else
    if(SL>0d0)then
        !icase=1
        ro=rl
        uo=ul
        Ptoto=Ptotl
        etoto=etotl
    else if(ustar>0d0)then
        !icase=2
        ro=rstarl
        uo=ustar
        Ptoto=Ptotstar
        etoto=etotstarl
    else if (SR>0d0)then
        !icase=3
        ro=rstarr
        uo=ustar
        Ptoto=Ptotstar
        etoto=etotstarr
    else
        !icase=4
        ro=rr
        uo=ur
        Ptoto=Ptotr
        etoto=etotr
    end if
#endif
    ! Compute the Godunov flux
    flux(i,r1,j1) = ro*uo
    flux(i,e1,j1) = (etoto+Ptoto)*uo
    flux(i,v1,j1) = ro*uo*uo+Ptoto
    flux(i,v2,j1) = ro*uo*merge(vl,vr,uo>0)
    flux(i,v3,j1) = ro*uo*merge(wl,wr,uo>0)
  end do
END SUBROUTINE hllc

!===============================================================================
!> Use the divergence of the HD fluxes to update the conservative variables.
!===============================================================================
SUBROUTINE div_flux (self, l, u, uold, unew, flux, dtds)
#ifdef OFFLOAD
  !$omp declare target
#endif
  class(mhd_t):: self
  integer:: l(4), u(4)
  real, contiguous:: uold(:,:,:,:), unew(:,:,:,:), flux(:,:,:,:,:)
  real:: dtds
  !.............................................................................
  real:: f, dmin0, dmin1, g0, g1, ek0, ek1, eth0, eth1, eth_min
  real:: ux, uy, uz, d0, d1, pg0, pg1, dtds1
  integer:: i1, i2, i3, i4, n_smalld, n_smallp
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%div_flux', itimer=itimer)
#ifdef OFFLOAD
  !$omp parallel do collapse(3) schedule(static,1)
#endif
  g0 = gamma
  g1 = merge(1d0,gamma-1d0,gamma==1d0)
  dmin0 = minval(uold(l(1):u(1),l(2):u(2),l(3):u(3),1))
  dtds1 = dtds
  call divfl
  do while (f < 1.0-courant_d)
    self%dtime = f*self%dtime
    dtds1      = f*dtds1
    if (verbose > 0) then
      write (stdlog,'(a,i6,f7.3,1p,3(2x,2e12.3))') &
        'div_flux: id,f,d0,d1,pg0,pg1,dtds_f=', &
        self%id, f, d0, d1, pg0, pg1, dtds1/dtds
    end if
    call divfl
  end do
  !.............................................................................
  if (verbose > 0) then
    unit = merge (stdout, stdlog, verbose > 1)
    if (n_smalld > 0) &
      write (unit,'(a,2i6,1p,g14.5)') &
        'mhd_t%div_flux: id,n_smalld,time=', self%id, n_smalld, self%time
    if (n_smallp > 0) &
      write (unit,'(a,2i6,1p,g14.5)') &
        'mhd_t%div_flux: id,n_smallp,time=', self%id, n_smallp, self%time
  end if
  !.............................................................................
  if (detailed_timer) call trace%end (itimer)
contains
!===============================================================================
!===============================================================================
subroutine divfl
  n_smalld = 0
  n_smallp = 0
  f = 1.0
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i4=1,5
        !$omp simd
        do i1=l(1),u(1)
          unew(i1,i2,i3,i4) = &
          uold(i1,i2,i3,i4) + dtds1*( &
            flux(i1,i2,i3,i4,1) - flux(i1+1,i2  ,i3  ,i4,1) &
          + flux(i1,i2,i3,i4,2) - flux(i1  ,i2+1,i3  ,i4,2) &
          + flux(i1,i2,i3,i4,3) - flux(i1  ,i2  ,i3+1,i4,3))
        end do
      end do
      do i1=l(1),u(1)
        d0 = uold(i1,i2,i3,1)
        ux = uold(i1,i2,i3,2)/d0
        uy = uold(i1,i2,i3,3)/d0
        uz = uold(i1,i2,i3,4)/d0
        ek0 = 0.5*d0*(ux**2+uy**2+uz**2)
        eth0 = uold(i1,i2,i3,5)
        d1 = unew(i1,i2,i3,1)
        ux = unew(i1,i2,i3,2)/d1
        uy = unew(i1,i2,i3,3)/d1
        uz = unew(i1,i2,i3,4)/d1
        n_smalld = merge (n_smalld+1, n_smalld, d1<smallr)
        d1 = max(d1,smallr)
        ek1 = 0.5*d1*(ux**2+uy**2+uz**2)
        eth_min = smallc**2/(g1*g0)*d1**g0
        eth1 = unew(i1,i2,i3,5) + (ek0-ek1)
        eth_min = smallc**2/(g1*g0)*d1**g0
        n_smallp = merge (n_smallp+1, n_smallp, eth1<eth_min)
        eth1 = max(eth1,eth_min)
        unew(i1,i2,i3,1) = d1
        unew(i1,i2,i3,5) = eth1
        f = min(f,d1/d0,eth1/eth0)
      end do
    end do
  end do
  f = exp(-(1.0-f))
end subroutine divfl
END SUBROUTINE div_flux

!===============================================================================
!===============================================================================
SUBROUTINE sources (l, u, unew, dt)
  integer :: l(3), u(3)
  real    :: unew(:,:,:,:)
  real    :: dt
  !.............................................................................
  integer :: i1, i2, i3
  real    :: d
  integer :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%forces', itimer=itimer)
  if (associated(scratch%force_per_unit_volume)) then
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
    do i3=l(3),u(3)
    do i2=l(2),u(2)
      !$omp simd
      do i1=l(1),u(1)
        !-----------------------------------------------------------------------
        ! d_ekin keeps track of the change of kinetic energy caused by the
        ! force, and adjust the total energy correspondingly
        !-----------------------------------------------------------------------
        unew(i1,i2,i3,u1) = &
        unew(i1,i2,i3,u1) + scratch%force_per_unit_volume(i1,i2,i3,1)*dt
        unew(i1,i2,i3,u2) = &
        unew(i1,i2,i3,u2) + scratch%force_per_unit_volume(i1,i2,i3,2)*dt
        unew(i1,i2,i3,u3) = &
        unew(i1,i2,i3,u3) + scratch%force_per_unit_volume(i1,i2,i3,3)*dt
      end do
    end do
    end do
  end if
  if (associated(scratch%force_per_unit_mass)) then
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
    do i3=l(3),u(3)
    do i2=l(2),u(2)
      !$omp simd
      do i1=l(1),u(1)
        d = unew(i1,i2,i3,r1)*dt
        unew(i1,i2,i3,u1) = &
        unew(i1,i2,i3,u1) + scratch%force_per_unit_mass(i1,i2,i3,1)*d
        unew(i1,i2,i3,u2) = &
        unew(i1,i2,i3,u2) + scratch%force_per_unit_mass(i1,i2,i3,2)*d
        unew(i1,i2,i3,u3) = &
        unew(i1,i2,i3,u3) + scratch%force_per_unit_mass(i1,i2,i3,3)*d
      end do
    end do
    end do
  end if
  if (associated(scratch%heating_per_unit_volume)) then
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
    do i3=l(3),u(3)
    do i2=l(2),u(2)
      !$omp simd
      do i1=l(1),u(1)
        !-----------------------------------------------------------------------
        ! d_ekin keeps track of the change of kinetic energy caused by the
        ! force, and adjust the total energy correspondingly
        !-----------------------------------------------------------------------
        unew(i1,i2,i3,e1) = &
        unew(i1,i2,i3,e1) + scratch%heating_per_unit_volume(i1,i2,i3)*dt
      end do
    end do
    end do
  end if
  if (associated(scratch%heating_per_unit_mass)) then
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
    do i3=l(3),u(3)
    do i2=l(2),u(2)
      !$omp simd
      do i1=l(1),u(1)
        d = unew(i1,i2,i3,r1)*dt
        unew(i1,i2,i3,e1) = &
        unew(i1,i2,i3,e1) + scratch%heating_per_unit_mass(i1,i2,i3)*d
      end do
    end do
    end do
  end if
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE sources

END MODULE mhd_mod
