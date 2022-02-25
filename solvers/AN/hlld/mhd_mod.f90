!===============================================================================
!> HLLD solver for GPUs (and CPUs).  To avoid potential issues with derived
!> data types, modules are (for now) only used to hold data and procedures
!>
!> call sequence inside update():
!>   ctoprim     :  compute primitive variables
!>   slopes      :  compute slopes of primitive variables
!>   predict     :  compute time derivatives and half-step primitives
!>   solve_1d    :  compute face values and solve for HD fluxes
!>   solve_2d    :  compute corner values and solver for EMFs
!>   div_flux    :  update conserved variables with fluxes and EMFs
!>
!> To minimize memory motions (also internally on the GPU), the left and right
!> face value arrays are 2-dimensional, with running index over the 1st dim,
!> which allows the 1-D solver to parallelize over the two redundant dimensions,
!> based only on access to the primitive variables and their (half) slopes.
!>
!> The results is walking through memory 7 times:
!> 1. ctoprim     uold -> prim
!> 2. slopes      prim -> dprm
!> 3. predict     prim -> prim, emf
!> 4. solve_1d    prim,dprim -> flux
!> 5. solve_2d    prim,emf   -> emf
!> 6. divflux     flux,uold  -> unew
!> 7. curlemf     emf ,uold  -> unew
!===============================================================================
MODULE mhd_mod
  USE io_mod
  USE trace_mod
  USE omp_lock_mod
  USE extras_mod
  USE eos_mod
  USE index_mod
  USE scratch_mod
  USE divb_clean_mod
  implicit none
  private
  type, public, extends(extras_t):: mhd_t
    real, pointer, contiguous, dimension(:,:,:,:)   :: prim, emf, bf
    real, pointer, contiguous, dimension(:,:,:,:,:) :: dprim, flux
    real, pointer, contiguous, dimension(:,:,:)     :: face, edge
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
    procedure:: post_update
    procedure:: update
    procedure, nopass:: ctoprim
  end type
  integer:: r1=1, u1=2, u2=3, u3=4, e1=5, m1=6, m2=7, m3=8
  real, parameter:: one=1d0, half=0.5, zero=0.0
#ifdef OFFLOAD
  !$omp declare target(n_ctoprim,n_solve_1d,n_solve_2d,n_predict,n_slopes)
  !$omp declare target(slope_type)
#else
  !save:: prim, gpu1, gpu2, gpu3, gpu4, gpu5, emf
  !!omp threadprivate (prim,gpu1,gpu2,gpu3,gpu5,emf)
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
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE pre_init (self)
  class(mhd_t), target:: self
  integer:: iostat
  namelist /an_mhd_params/ verbose, courant, courant_type, slope_type, csound, &
                           courant_d, smallr, smallc, unsigned, detailed_timer
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
    write (io%namelist,an_mhd_params)
    gamma = eos%gamma
  end if
  call io%lock%unset
  self%courant = courant
  self%ng = 3
  self%nt = 5
  self%nv = 8
  self%nw = 1
  self%kind  = 'AN_hlld'
  self%type  = 'mhd_t'
  self%etype = 'thermal'
  self%staggered = .false.
  call self%idx%init (energy=5, mhd=.true.)
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
    self%mesh(i)%h(self%idx%bx+i-1) = -0.5
  end do
  self%unsigned(self%idx%d) = unsigned
  self%unsigned(self%idx%e) = unsigned
  self%pervolume(self%idx%px) = unsigned
  self%pervolume(self%idx%py) = unsigned
  self%pervolume(self%idx%pz) = unsigned
  call divb%init
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
  m = shape(self%mem)
  m(4) = self%nv
  allocate (self%dprim(m(1),m(4),3,m(2),m(3)))
  allocate (self%flux (m(1),m(4),3,m(2),m(3)))
  allocate (self%emf  (m(1),m(2),  m(3),3))
  allocate (self%prim (m(1),m(4),  m(2),m(3)))
  allocate (self%face (m(1),m(4),2))
  allocate (self%edge (m(1),m(4),4))
  !-----------------------------------------------------------------------------
  uold => self%mem(:,:,:,:,self%it,1)
  unew => self%mem(:,:,:,:,self%new,1)
  self%bf => uold(:,:,:,self%idx%bx:self%idx%bz)
  l(4) = 1
  u(4) = self%nv
  !-----------------------------------------------------------------------------
  ! The primitive values are needed in one cell more than the slopes, and to
  ! compute them requires access to one more cell below, so from li-3
  !-----------------------------------------------------------------------------
  l(1:3) = self%mesh%lb
  u(1:3) = self%mesh%ub
  call ctoprim (self, l, u, uold, self%prim, u_max)
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
!> In the sketch of guard zones below, 1-3 are the guard cells, + are internal
!> 1---2---3---+................+---3---2---1        ! full range
!>   2---3---+................+---3---2---1          ! B
!>     2---3---+................+---3---2            ! ctoprim valid
!>         3---+................+---3                ! slopes valid
!>           3---+................+---3              ! right valid
!>       3---+................+---3                  ! left valid
!>           3---+................+---3              ! lft valid
!>           3---+................+                  ! riemann valid
!>             +................+                    ! div(flux) valid
!===============================================================================
SUBROUTINE update (self)
#ifdef OFFLOAD
  !$omp declare target
#endif
  class(mhd_t)              :: self
  integer                   :: l(4), u(4), j
  real, contiguous, pointer :: uold(:,:,:,:), unew(:,:,:,:)
  real,             pointer :: b(:,:,:,:)
  type(lock_t),     pointer :: lock
  real                      :: dtds, dt2, dt
  integer,             save :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Use one of several methods to clean div(B) from guard zones
  !-----------------------------------------------------------------------------
  b => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,1)
  call divb%clean (self%mesh, b)
  !-----------------------------------------------------------------------------
  ! Possibly do output
  !-----------------------------------------------------------------------------
  call self%output
  !-----------------------------------------------------------------------------
  ! With these l:u bounds, the slopes are good in the whole l:u domaini, except
  ! the B-component slopes along the component direction should not be used at
  ! the uppermost index (u = ub-1)
  !-----------------------------------------------------------------------------
  l(1:3) = self%mesh%lb+1
  u(1:3) = self%mesh%ub-1
  l(4) = 1
  u(4) = self%nv
  call slopes   (l, u, self%prim, self%dprim)
  !-----------------------------------------------------------------------------
  ! The predictor needs to compute values one cell lower than the flux, to allow
  ! solve_1d and solve_2d to reference cells displaced by -1
  !-----------------------------------------------------------------------------
  l(1:3) = self%mesh%li-1
  u(1:3) = self%mesh%uo
  dtds = self%dtime/self%ds(1)
  dt2 = self%dtime/2.0
  call predict (l, u, self%prim, self%dprim, self%bf, self%emf, &
                dtds, dt2)
  !-----------------------------------------------------------------------------
  ! The fluxes need to be computed in the interval li-uo, to allow dd_up
  !-----------------------------------------------------------------------------
  l(1:3) = self%mesh%li-1
  u(1:3) = self%mesh%uo
  call solve_1d (l, u, self%prim, self%dprim, self%flux, self%face)
  call solve_2d (l, u, self%prim, self%dprim, self%emf, self%edge)
  !.............................................................................
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
  call div_flux (self, l, u, uold, unew, self%flux, self%emf, dtds)
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
  deallocate (self%prim)
  deallocate (self%dprim)
  deallocate (self%flux)
  deallocate (self%emf)
  call self%extras_t%post_update ()
  call self%gpatch_t%counter_update ()
  !.............................................................................
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE post_update

!===============================================================================
!> Conservative to primitive variables. Like all procedures below, it does not
!> refer to any structure components -- information is passed via argument lists
!===============================================================================
SUBROUTINE ctoprim (self, l, u, uold, prim, u_max)
#ifdef OFFLOAD
  !$omp declare target
#endif
  class(mhd_t)     :: self
  integer          :: l(4), u(4)
  real, contiguous :: uold(:,:,:,:), prim(:,:,:,:)
  real             :: u_max, u_fast, u_gas_1, u_gas_2, cs2 
  !.............................................................................
  integer:: i1, i2, i3, n_smalld, n_smallp, unit
  real:: smallp, g0, g1
  real:: d, e, px, py, pz, ux, uy, uz, bx, by, bz, pg, u2_max, b2, pg0, pg1
  integer, save:: itimer=0
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
  cs2 = csound**2
  !-----------------------------------------------------------------------------
  ! This loop leaves OK values in the whole domain, except that each of the
  ! B-components has a fill-in value at the upper boundary.
  !-----------------------------------------------------------------------------
  n_smalld = 0
  n_smallp = 0
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !$omp simd
    do i1=l(1),u(1)
      d  = uold(i1,i2,i3,r1)
      e  = uold(i1,i2,i3,e1)
      px = uold(i1,i2,i3,u1)
      py = uold(i1,i2,i3,u2)
      pz = uold(i1,i2,i3,u3)
      n_smalld = merge(n_smalld+1,n_smalld,d<smallr)
      d = max(d,smallr)
      prim(i1,r1,i2,i3) = log(d)
      ux = px/d
      uy = py/d
      uz = pz/d
      prim(i1,u1,i2,i3) = ux
      prim(i1,u2,i2,i3) = uy
      prim(i1,u3,i2,i3) = uz
      bx = (uold(    i1        ,i2,i3,m1)  &
         +  uold(min(i1+1,u(1)),i2,i3,m1))*0.5
      by = (uold(i1,    i2        ,i3,m2)  &
         +  uold(i1,min(i2+1,u(2)),i3,m2))*0.5
      bz = (uold(i1,i2,    i3        ,m3)  &
         +  uold(i1,i2,min(i3+1,u(3)),m3))*0.5
      prim(i1,m1,i2,i3) = bx
      prim(i1,m2,i2,i3) = by
      prim(i1,m3,i2,i3) = bz
      b2 = bx**2+by**2+bz**2
      pg0 = g1*e
      pg1 = cs2*d
      pg = merge (pg1, pg0, gamma==1d0)
      !--- this scales correctly with density^gamma, but assumes average d=1 ---
      smallp = smallc**2/g0*d**g0
      n_smallp = merge(n_smallp+1,n_smallp,pg<smallp)
      pg = max(pg,smallp)
      u_fast = sqrt((gamma*pg+b2)/d)
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
        'mhd_t%ctoprim: id,n_smalld,time=', self%id, n_smalld, self%time
    if (n_smallp > 0) &
      write (unit,'(a,2i6,1p,g14.5)') &
        'mhd_t%ctoprim: id,n_smallp,time=', self%id, n_smallp, self%time
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
  integer, save:: itimer=0
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
SUBROUTINE predict (l, u, q, dq, bf, e, dtds, dt2)
#ifdef OFFLOAD
  !$omp declare target
#endif
  integer                    :: l(4), u(4)
  real, dimension(:,:,:,:)   :: q      ! cell center
  real, dimension(:,:,:,:,:) :: dq     ! half slopes
  real, dimension(:,:,:,:)   :: bf     ! cell faces
  real, dimension(:,:,:,:)   :: e      ! electric field
  real                       :: dtds   ! dt/ds
  real                       :: dt2    ! dt/2
  intent(in)                 :: dq
  contiguous                 :: q, dq, e
  !.............................................................................
  integer :: i1, i2, i3, i4
  real    :: g0, c, dtds2
  real    :: ld, d_ld, ld_x, ld_y, ld_z, div_u
  real    :: lp, d_lp, lp_x, lp_y, lp_z, tt
  real    :: ux, ux_x, ux_y, ux_z
  real    :: uy, uy_x, uy_y, uy_z
  real    :: uz, uz_x, uz_y, uz_z
  real    :: bx,       bx_y, bx_z, ex
  real    :: by, by_x,       by_z, ey
  real    :: bz, bz_x, bz_y,       ez
  real    :: ff(l(1):u(1),3), qq(l(1):u(1))
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%predict', itimer=itimer)
  g0 = gamma
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
  dtds2 = dtds*0.5
  e = 0.0
  do i3=l(3),u(3)+2
  do i2=l(2),u(2)+2
    !$omp simd
    do i1=l(1),u(1)+2
      !-------------------------------------------------------------------------
      ! Down-down edge centrered electric field at initial time.  The max range
      ! of this loop is lb+1:ub, because it does not rely on slopes of U.
      !-------------------------------------------------------------------------
      ux   =  q(i1,u1,i2,i3)
      uy   =  q(i1,u2,i2,i3)
      uz   =  q(i1,u3,i2,i3)
      ux_x = (q(i1,u1,i2,i3)-q(i1-1,u1,i2  ,i3  ))*0.5
      ux_y = (q(i1,u1,i2,i3)-q(i1  ,u1,i2-1,i3  ))*0.5
      ux_z = (q(i1,u1,i2,i3)-q(i1  ,u1,i2  ,i3-1))*0.5
      uy_x = (q(i1,u2,i2,i3)-q(i1-1,u2,i2  ,i3  ))*0.5
      uy_y = (q(i1,u2,i2,i3)-q(i1  ,u2,i2-1,i3  ))*0.5
      uy_z = (q(i1,u2,i2,i3)-q(i1  ,u2,i2  ,i3-1))*0.5
      uz_x = (q(i1,u3,i2,i3)-q(i1-1,u3,i2  ,i3  ))*0.5
      uz_y = (q(i1,u3,i2,i3)-q(i1  ,u3,i2-1,i3  ))*0.5
      uz_z = (q(i1,u3,i2,i3)-q(i1  ,u3,i2  ,i3-1))*0.5
      !-------------------------------------------------------------------------
      ! The magnetic field components here are staggered in their own direction
      ! and are only used to compute the EMF = - (u x B)
      !-------------------------------------------------------------------------
      bx   =  bf(i1,i2,i3,1)
      by   =  bf(i1,i2,i3,2)
      bz   =  bf(i1,i2,i3,3)
      !-------------------------------------------------------------------------
      ! These slopes, perpendicular to the component direction, return values
      ! that are edge-centered (down) in the component and slope directions
      !-------------------------------------------------------------------------
      bx_y = (bf(i1,i2,i3,1)-bf(i1,i2-1,i3,1))*0.5
      bx_z = (bf(i1,i2,i3,1)-bf(i1,i2,i3-1,1))*0.5
      by_z = (bf(i1,i2,i3,2)-bf(i1,i2,i3-1,2))*0.5
      by_x = (bf(i1,i2,i3,2)-bf(i1-1,i2,i3,2))*0.5
      bz_x = (bf(i1,i2,i3,3)-bf(i1-1,i2,i3,3))*0.5
      bz_y = (bf(i1,i2,i3,3)-bf(i1,i2-1,i3,3))*0.5
      ex = -(uy-uy_y-uy_z)*(bz-bz_y) + (uz-uz_z-uz_y)*(by-by_z)
      ey = -(uz-uz_z-uz_x)*(bx-bx_z) + (ux-ux_x-ux_z)*(bz-bz_x)
      ez = -(ux-ux_x-ux_y)*(by-by_x) + (uy-uy_y-uy_x)*(bx-bx_y)
      e(i1,i2,i3,1) = dtds2*ex
      e(i1,i2,i3,2) = dtds2*ey
      e(i1,i2,i3,3) = dtds2*ez
    end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Include external source terms in the predictor
  !-----------------------------------------------------------------------------
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
  do i3=l(3),u(3)
  do i2=l(2),u(2)
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
      ld = ld - dtds*(ux*ld_x + uy*ld_y + uz*ld_z +    div_u)
      lp = lp - dtds*(ux*lp_x + uy*lp_y + uz*lp_z + g0*div_u) + qq(i1)
      !-------------------------------------------------------------------------
      ! Velocity evolution
      !-------------------------------------------------------------------------
      tt = exp(lp-ld)
      ux = ux - dtds*(ux*ux_x + uy*ux_y + uz*ux_z + tt*lp_x) + ff(i1,1)
      uy = uy - dtds*(ux*uy_x + uy*uy_y + uz*uy_z + tt*lp_y) + ff(i1,2)
      uz = uz - dtds*(ux*uz_x + uy*uz_y + uz*uz_z + tt*lp_z) + ff(i1,3)
      !-------------------------------------------------------------------------
      ! Store results
      !-------------------------------------------------------------------------
      q(i1,r1,i2,i3) = ld
      q(i1,e1,i2,i3) = lp
      q(i1,u1,i2,i3) = ux
      q(i1,u2,i2,i3) = uy
      q(i1,u3,i2,i3) = uz
    end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Predict the B-field at t+dt/2, for indices l:u+1, with input required in the
  ! index range l:u+2.  The max index range for these loops is lb+1:ub-2, since
  ! the E-field is good to ub-1.
  !-----------------------------------------------------------------------------
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
  do i3=l(3),u(3)+1
  do i2=l(2),u(2)+1
    !$omp simd
    do i1=l(1),u(1)+1
      !-------------------------------------------------------------------------
      ! Magnetic field evolution; dB/dt = -curl(E), this evolves the down-
      ! staggered magnetic field components at cell interfaces
      !-------------------------------------------------------------------------
      q(i1,m1,i2,i3) = &
      bf(i1,i2,i3,1) - (e(i1,i2+1,i3,3) - e(i1,i2,i3,3)) &
                     + (e(i1,i2,i3+1,2) - e(i1,i2,i3,2))
      q(i1,m2,i2,i3) = &
      bf(i1,i2,i3,2) - (e(i1,i2,i3+1,1) - e(i1,i2,i3,1)) &
                     + (e(i1+1,i2,i3,3) - e(i1,i2,i3,3))
      q(i1,m3,i2,i3) = &
      bf(i1,i2,i3,3) - (e(i1+1,i2,i3,2) - e(i1,i2,i3,2)) &
                     + (e(i1,i2+1,i3,1) - e(i1,i2,i3,1))
    end do
  end do
  end do
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !$omp simd
    do i1=l(1),u(1)
      !-------------------------------------------------------------------------
      ! Compute predicted cell-centered B-components; these are good in the
      ! index range up to uo = ub-2, except for the coordinate direction, wich
      ! only has good values up to ui = ub-3
      !-------------------------------------------------------------------------
      q(i1,m1,i2,i3) = (q(i1,m1,i2,i3)+q(i1+1,m1,i2,i3))*0.5
      q(i1,m2,i2,i3) = (q(i1,m2,i2,i3)+q(i1,m2,i2+1,i3))*0.5
      q(i1,m3,i2,i3) = (q(i1,m3,i2,i3)+q(i1,m3,i2,i3+1))*0.5
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
!>
!> Two 2-D arrays are passed to the actual Riemann solver, with the x-index in
!> the 1st slot and the variable index in the 2nd.  The Riemann solver needs to
!> pick up the parallel and perpendicular directions, based on the direction
!> index j, passed as an argument.
!===============================================================================
SUBROUTINE solve_1d (l, u, q, dq, flux, a)
#ifdef OFFLOAD
  !$omp declare target
#endif
  integer          :: l(4), u(4)
  real, contiguous :: q(:,:,:,:), dq(:,:,:,:,:), flux(:,:,:,:,:), a(:,:,:)
  intent(in)       :: q, dq
  !.............................................................................
  integer :: i1, i2, i3, i4, j, k1, k2, k3, kk(3), iv
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%solve_1d', itimer=itimer)
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    do j=1,3
      kk(:) = 0
      kk(j) = 1
      k1 = kk(1)
      k2 = kk(2)
      k3 = kk(3)
      do i4=l(4),u(4)
        !$omp simd
        do i1=l(1),u(1)
          a(i1,i4,1) =  q(i1-k1,i4,  i2-k2,i3-k3)  & ! qleft
                     + dq(i1-k1,i4,j,i2-k2,i3-k3)    ! ds/2 increment
          a(i1,i4,2) =  q(i1   ,i4,  i2   ,i3   )  & ! qright
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
      call hlld_1d (l, u, a(:,:,1), a(:,:,2), flux(:,:,:,i2,i3), j)
    end do
  end do
  end do
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE solve_1d

!===============================================================================
!> The RAMSES 1-D HLLD solver, verbatim.
!===============================================================================
SUBROUTINE hlld_1d (lb, ub, qleft, qright, flx, j1)
#ifdef OFFLOAD
  !$omp declare target
#endif
  integer:: lb(4), ub(4)
  real, dimension(:,:), intent(inout) :: qleft
  real, dimension(:,:), intent(in   ) :: qright
  real, dimension(:,:,:), intent(out) :: flx
  integer:: j1
  !.............................................................................
  real::SL,SR
  real::entho,lim
  real::rl,pl,ul,ecinl,etotl,ptotl,Bl,Cl,wstarl,vtotBl,emagl,SAL,sqrrstarl,vl,wl
  real::rr,pr,ur,ecinr,etotr,ptotr,Br,Cr,wstarr,vtotBr,emagr,SAR,sqrrstarr,vr,wr
  real::cfastl,rcl,rstarl,el,vdotBstarl,vstarl,etotstarstarl,etotstarl,vdotBl
  real::cfastr,rcr,rstarr,er,vdotBstarr,vstarr,etotstarstarr,etotstarr,vdotBr
  real::ustar,ptotstar,estar,vstarstar,wstarstar,ro,uo,ptoto,etoto,vo,wo
  real::A,sgnm
  real::Bstarl,Bstarr,Cstarl,Cstarr,calfvenl,calfvenr
  real::Bo,Co,Bstarstar,Cstarstar,vdotBstarstar,vdotBo,emago
  real::d,P,u,v,w,B,C,B2,c2,d2,cf
  logical:: flag, illconditioned
  integer::ivar,i,icase, ndiff, j2, j3, v1, v2, v3, n1, n2, n3
  !-----------------------------------------------------------------------------
  entho = merge(1d0, 1d0/max(gamma-1d0,1d-6), gamma==1d0)
  j2 = 1 + mod(j1  ,3)
  j3 = 1 + mod(j1+1,3)
  v1 = u1+j1-1
  v2 = u1+j2-1
  v3 = u1+j3-1
  n1 = m1+j1-1
  n2 = m1+j2-1
  n3 = m1+j3-1
!define PRIVATE
#ifdef PRIVATE
  !$omp simd private(sgnm,d,P,u,v,w,A,B,C,B2,c2,d2,cf,Ptotstar,estar, &
  !$omp &  rl,vl,ecinl,emagl,etotl,Ptotl,vdotBl,Bl,Cl,cfastl,SL, &
  !$omp &  rr,vr,ecinr,emagr,etotr,Ptotr,vdotBr,Br,Cr,cfastr,SR, &
  !$omp &  rcl,rstarl,el,Bstarl,Cstarl,vstarl,wstarl, &
  !$omp &  rcr,rstarr,er,Bstarr,Cstarr,vstarr,wstarr, &
  !$omp &  vdotBstarl,etotstarl,sqrrstarl,calfvenl,SAR, &
  !$omp &  vdotBstarr,etotstarr,sqrrstarr,calfvenr,SAL, &
  !$omp &  vstarstar,wstarstar,Bstarstar,Cstarstar,vdotBstarstar, &
  !$omp &  icase,ro,uo,vo,wo,Bo,Co,Ptoto,etoto,vdotBo,emago,ustar)
#else
  !$omp simd
#endif
  do i=lb(1),ub(1)
    !
    ! Enforce continuity of normal component
    A=half*(qleft(i,n1)+qright(i,n1))
    sgnm=sign(one,A)
    !
    ! Left variables
    rl=qleft(i,r1)
    Pl=qleft(i,e1)
    ul=qleft(i,v1); vl=qleft(i,v2); wl=qleft(i,v3)
                    Bl=qleft(i,n2); Cl=qleft(i,n3)
    ecinl = half*(ul*ul+vl*vl+wl*wl)*rl
    emagl = half*(A*A+Bl*Bl+Cl*Cl)
    etotl = Pl*entho+ecinl+emagl
    Ptotl = Pl + emagl
    vdotBl= ul*A+vl*Bl+wl*cl
    !
    ! Right variables
    rr=qright(i,r1)
    Pr=qright(i,e1)
    ur=qright(i,v1); vr=qright(i,v2); wr=qright(i,v3)
                     Br=qright(i,n2); Cr=qright(i,n3)
    ecinr = half*(ur*ur+vr*vr+wr*wr)*rr
    emagr = half*(A*A+Br*Br+Cr*Cr)
    etotr = Pr*entho+ecinr+emagr
    Ptotr = Pr + emagr
    vdotBr= ur*A+vr*Br+wr*Cr
    !
    ! Find the largest eigenvalues in the normal direction to the interface
    B2 = A*A+Bl*Bl+Cl*Cl
    c2 = gamma*Pl/rl
    d2 = half*(B2/rl+c2)
    cfastl = sqrt(d2 + sqrt(max(d2*d2-c2*A*A/rl,zero)) )
    B2 = A*A+Br*Br+Cr*Cr
    c2 = gamma*Pr/rr
    d2 = half*(B2/rr+c2)
    cfastr = sqrt(d2 + sqrt(max(d2*d2-c2*A*A/rr,zero)) )
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
    estar =rl*(SL-ul)*(SL-ustar)-A**2
    el    =rl*(SL-ul)*(SL-ul   )-A**2
    !
    if(abs(estar) <= 1e-4*A**2)then
      vstarl=vl
      Bstarl=Bl
      wstarl=wl
      Cstarl=Cl
    else
      vstarl=vl-A*Bl*(ustar-ul)/estar
      Bstarl=Bl*el/estar
      wstarl=wl-A*Cl*(ustar-ul)/estar
      Cstarl=Cl*el/estar
    endif
    vdotBstarl=ustar*A+vstarl*Bstarl+wstarl*Cstarl
    etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar+A*(vdotBl-vdotBstarl))/(SL-ustar)
    sqrrstarl=sqrt(rstarl)
    calfvenl=abs(A)/sqrrstarl
    SAL=ustar-calfvenl
    !
    ! Right star region variables
    rstarr=rr*(SR-ur)/(SR-ustar)
    estar =rr*(SR-ur)*(SR-ustar)-A**2
    er    =rr*(SR-ur)*(SR-ur   )-A**2
#ifdef MERGE
    illconditioned = abs(estar) < 1d-4*A**2
    vstarr = merge (vr, vr-A*Br*(ustar-ur)/estar, illconditioned)
    Bstarr = merge (Br, Br*er/estar             , illconditioned)
    wstarr = merge (wr, wr-A*Cr*(ustar-ur)/estar, illconditioned)
    cstarr = merge (cr, Cr*er/estar             , illconditioned)
#else
  if(abs(estar) <= 1e-4*A**2)then
      vstarr=vr
      Bstarr=Br
      wstarr=wr
      Cstarr=Cr
  else
      vstarr=vr-A*Br*(ustar-ur)/estar
      Bstarr=Br*er/estar
      wstarr=wr-A*Cr*(ustar-ur)/estar
      Cstarr=Cr*er/estar
  endif
#endif
    !
    vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr
    etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar+A*(vdotBr-vdotBstarr))/(SR-ustar)
    sqrrstarr=sqrt(rstarr)
    calfvenr=abs(A)/sqrrstarr
    SAR=ustar+calfvenr
    !
    ! Double star region variables
    vstarstar=(sqrrstarl*vstarl+sqrrstarr*vstarr+sgnm*(Bstarr-Bstarl))/(sqrrstarl+sqrrstarr)
    wstarstar=(sqrrstarl*wstarl+sqrrstarr*wstarr+sgnm*(Cstarr-Cstarl))/(sqrrstarl+sqrrstarr)
    Bstarstar=(sqrrstarl*Bstarr+sqrrstarr*Bstarl+sgnm*sqrrstarl*sqrrstarr*(vstarr-vstarl))/(sqrrstarl+sqrrstarr)
    Cstarstar=(sqrrstarl*Cstarr+sqrrstarr*Cstarl+sgnm*sqrrstarl*sqrrstarr*(wstarr-wstarl))/(sqrrstarl+sqrrstarr)
    vdotBstarstar=ustar*A+vstarstar*Bstarstar+wstarstar*Cstarstar
    etotstarstarl=etotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
    etotstarstarr=etotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)
    !
    ! Sample the solution at x/t=0
    if(SL > 0.0)then
      icase=1
      ro=rl
      uo=ul
      vo=vl
      wo=wl
      Bo=Bl
      Co=Cl
      Ptoto=Ptotl
      etoto=etotl
      vdotBo=vdotBl
      emago=emagl
    else if(SAL > 0.0)then
      icase=2
      ro=rstarl
      uo=ustar
      vo=vstarl
      wo=wstarl
      Bo=Bstarl
      Co=Cstarl
      Ptoto=Ptotstar
      etoto=etotstarl
      vdotBo=vdotBstarl
      emago=half*(A*A+Bo*Bo+Co*Co)
    else if(ustar > 0.0)then
      icase=3
      ro=rstarl
      uo=ustar
      vo=vstarstar
      wo=wstarstar
      Bo=Bstarstar
      Co=Cstarstar
      Ptoto=Ptotstar
      etoto=etotstarstarl
      vdotBo=vdotBstarstar
      emago=half*(A*A+Bo*Bo+Co*Co)
    else if(SAR > 0.0)then
      icase=4
      ro=rstarr
      uo=ustar
      vo=vstarstar
      wo=wstarstar
      Bo=Bstarstar
      Co=Cstarstar
      Ptoto=Ptotstar
      etoto=etotstarstarr
      vdotBo=vdotBstarstar
      emago=half*(A*A+Bo*Bo+Co*Co)
    else if (SR > 0.0)then
      icase=5
      ro=rstarr
      uo=ustar
      vo=vstarr
      wo=wstarr
      Bo=Bstarr
      Co=Cstarr
      Ptoto=Ptotstar
      etoto=etotstarr
      vdotBo=vdotBstarr
      emago=half*(A*A+Bo*Bo+Co*Co)
    else
      icase=6
      ro=rr
      uo=ur
      vo=vr
      wo=wr
      Bo=Br
      Co=Cr
      Ptoto=Ptotr
      etoto=etotr
      vdotBo=vdotBr
      emago=emagr
    end if
    !---------------------------------------------------------------------------
    ! Add a diffusive flux in the down grad(ln(rho)) direction, which is only
    ! nonzero when the smallest of the left/right density is within a factor 10
    ! from smallr, and which has at most magnitude of 0.25 times the largest of
    ! the left/right fast mode speeds.
    !---------------------------------------------------------------------------
#ifdef DIFF
    ustar = log(rl/rr)*max(0.0, log((10.*smallr)/min(rl,rr)))
    ustar = 0.25*max(cfastr,cfastl)*min(1.0, max(-1.0, ustar))
    uo = uo + ustar
#endif
    !
#ifdef DBG
    if (rprint > 0 .and. i > 11 .and. i < 21) then
      if (i == 12) print *, &
        '  i cse   ro      rcl     rcr     SL      SR      ul      ur'// &
        '      ustar   uo      Pl      Pr      Ptotl   Ptotr   Ptoto'
      print '(2i4,16f8.3)', &
        i,icase,ro,rcl,rcr,SL,SR,ul,ur,ustar,uo,Pl,Pr,Ptotl,Ptotr,Ptoto
      if (i == 20) rprint = 0
    end if
#endif
    ! Compute the Godunov fluxes
    flx(i,r1,j1) = ro*uo
    flx(i,e1,j1) = (etoto+Ptoto)*uo-A*vdotBo
    flx(i,v1,j1) = ro*uo*uo+Ptoto-A*A
    flx(i,v2,j1) = ro*uo*vo-A*Bo
    flx(i,v3,j1) = ro*uo*wo-A*Co
  end do
END SUBROUTINE hlld_1d

!===============================================================================
!===============================================================================
subroutine hlld_ramses_mhd (lb, ub, qleft, qright, flux, j1)
  IMPLICIT none
  integer:: lb(4), ub(4)
  real, dimension(:,:), intent(inout) :: qleft
  real, dimension(:,:), intent(in   ) :: qright
  real, dimension(:,:,:), intent(out) :: flux
  integer:: j1
  !.............................................................................
  real::sl,sr
  real::entho,lim
  real::rl,pl,ul,ecinl,etotl,ptotl,bl,cl,wstarl,vtotbl,emagl,sal,sqrrstarl,vl,wl
  real::rr,pr,ur,ecinr,etotr,ptotr,br,cr,wstarr,vtotbr,emagr,sar,sqrrstarr,vr,wr
  real::cfastl,rcl,rstarl,el,vdotbstarl,vstarl,etotstarstarl,etotstarl,vdotbl
  real::cfastr,rcr,rstarr,er,vdotbstarr,vstarr,etotstarstarr,etotstarr,vdotbr
  real::ustar,ptotstar,estar,vstarstar,wstarstar
  real::ro,uo,ptoto,etoto,vo,wo
  real::a,sgnm
  real::bstarl,bstarr,cstarl,cstarr,calfvenl,calfvenr
  real::bo,co,bstarstar,cstarstar,vdotbstarstar,vdotbo,emago
#define INLINE
#ifdef INLINE
  real::d, p, u, v, w, b, c, b2, c2, d2, cf
#endif
  integer ::ivar, i, icase, j2, j3, v1, v2, v3, n1, n2, n3
  logical:: flag, illconditioned
  integer, save:: itimer=0
  !.............................................................................
  !
  j2 = 1 + mod(j1,3)
  j3 = 1 + mod(j2,3)
  v1 = u1+j1-1
  v2 = u1+j2-1
  v3 = u1+j3-1
  n1 = u1+j1-1
  n2 = u1+j2-1
  n3 = u1+j3-1
  !
  ! constants
  if (gamma==1d0) then
    entho = one
  else
    entho = one/max(eos%gamma-one,1e-6)
  end if
  !
  do i=lb(1),ub(1)
    !
    ! Enforce continuity of normal component
    A=half*(qleft(i,n1)+qright(i,n1))
    sgnm=sign(one,A)
    !qleft(i,4)=A; qright(i,4)=A
    !
    ! Left variables
    rl=qleft(i,r1); Pl=qleft(i,e1); ul=qleft(i,v1)
    vl=qleft(i,v2); Bl=qleft(i,m2); wl=qleft(i,v3); Cl=qleft(i,m3)
    ecinl = half*(ul*ul+vl*vl+wl*wl)*rl
    emagl = half*(A*A+Bl*Bl+Cl*Cl)
    etotl = Pl*entho+ecinl+emagl
    Ptotl = Pl + emagl
    vdotBl= ul*A+vl*Bl+wl*cl
    !
    ! Right variables
    rr=qright(i,r1); Pr=qright(i,e1); ur=qright(i,v1)
    vr=qright(i,v2); Br=qright(i,m2); wr=qright(i,v3); Cr=qright(i,m3)
    ecinr = half*(ur*ur+vr*vr+wr*wr)*rr
    emagr = half*(A*A+Br*Br+Cr*Cr)
    etotr = Pr*entho+ecinr+emagr
    Ptotr = Pr + emagr
    vdotBr= ur*A+vr*Br+wr*Cr
    ! fgdnv(i,1) = ro*uo                                      ! rho
    ! fgdnv(i,2) = (etoto+Ptoto)*uo-A*vdotBo                  ! 1st U comp
    ! fgdnv(i,3) = ro*uo*uo+Ptoto+emago*(lim-1.0_dp)-A*A*lim  ! total E
    ! fgdnv(i,4) = zero                                       ! 1st B comp
    ! fgdnv(i,5) = ro*uo*vo-A*Bo*lim                          ! 2nd U comp
    ! fgdnv(i,6) = Bo*uo-A*vo                                 ! 2nd B comp
    ! fgdnv(i,7) = ro*uo*wo-A*Co*lim                          ! 3rd U comp
    ! fgdnv(i,8) = Co*uo-A*wo                                 ! 3d B comp
    !
    ! Find the largest eigenvalues in the normal direction to the interface
    d=qleft(i,r1); P=qleft(i,e1)
    u=qleft(i,v1); v=qleft(i,v2); w=qleft(i,v3)
    A=qleft(i,n1); B=qleft(i,n2); C=qleft(i,n3)
    B2 = A*A+B*B+C*C
    c2 = eos%gamma*P/d
    d2 = half*(B2/d+c2)
    cf = sqrt( d2 + sqrt(max(d2*d2-c2*A*A/d,zero)) )
    cfastl=cf
    !
    d=qright(i,r1); P=qright(i,e1)
    u=qright(i,v1); v=qright(i,v2); w=qright(i,v3)
    A=qright(i,n1); B=qright(i,n2); C=qright(i,n3)
    B2 = A*A+B*B+C*C
    c2 = eos%gamma*P/d
    d2 = half*(B2/d+c2)
    cf = sqrt( d2 + sqrt(max(d2*d2-c2*A*A/d,zero)) )
    cfastr=cf
    !
    ! Compute HLL wave speed
    SL=min(ul,ur)-max(cfastl,cfastr)
    SR=max(ul,ur)+max(cfastl,cfastr)
    !
    lim = 1d0
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
    estar =rl*(SL-ul)*(SL-ustar)-A**2
    el    =rl*(SL-ul)*(SL-ul   )-A**2
    !if (dbg) print 1,'hlld4', estar,rl,SL,ul,ustar,A
    if(abs(estar) < 1d-4*A**2)then
      vstarl=vl
      Bstarl=Bl
      wstarl=wl
      Cstarl=Cl
    else
      vstarl=vl-A*Bl*(ustar-ul)/estar
      Bstarl=Bl*el/estar
      wstarl=wl-A*Cl*(ustar-ul)/estar
      Cstarl=Cl*el/estar
    endif
     vdotBstarl=ustar*A+vstarl*Bstarl+wstarl*Cstarl
    etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar+A*(vdotBl-vdotBstarl))/(SL-ustar)
    sqrrstarl=sqrt(rstarl)
    calfvenl=abs(A)/sqrrstarl
    SAL=ustar-calfvenl
    !
    ! Right star region variables
    rstarr=rr*(SR-ur)/(SR-ustar)
    estar =rr*(SR-ur)*(SR-ustar)-A**2
    er    =rr*(SR-ur)*(SR-ur   )-A**2
  if(abs(estar) < 1d-4*A**2)then
      vstarr=vr
      Bstarr=Br
      wstarr=wr
      Cstarr=Cr
  else
      vstarr=vr-A*Br*(ustar-ur)/estar
      Bstarr=Br*er/estar
      wstarr=wr-A*Cr*(ustar-ur)/estar
      Cstarr=Cr*er/estar
  endif
    vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr
    etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar+A*(vdotBr-vdotBstarr))/(SR-ustar)
    sqrrstarr=sqrt(rstarr)
    calfvenr=abs(A)/sqrrstarr
    SAR=ustar+calfvenr
    !
    ! Double star region variables
    vstarstar=(sqrrstarl*vstarl+sqrrstarr*vstarr+sgnm*(Bstarr-Bstarl))/(sqrrstarl+sqrrstarr)
    wstarstar=(sqrrstarl*wstarl+sqrrstarr*wstarr+sgnm*(Cstarr-Cstarl))/(sqrrstarl+sqrrstarr)
    Bstarstar=(sqrrstarl*Bstarr+sqrrstarr*Bstarl+sgnm*sqrrstarl*sqrrstarr*(vstarr-vstarl))/(sqrrstarl+sqrrstarr)
    Cstarstar=(sqrrstarl*Cstarr+sqrrstarr*Cstarl+sgnm*sqrrstarl*sqrrstarr*(wstarr-wstarl))/(sqrrstarl+sqrrstarr)
    vdotBstarstar=ustar*A+vstarstar*Bstarstar+wstarstar*Cstarstar
    etotstarstarl=etotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
    etotstarstarr=etotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)
    !
    ! Sample the solution at x/t=0
    if(SL>0d0)then
       icase=1
       ro=rl
       uo=ul
       vo=vl
       wo=wl
       Bo=Bl
       Co=Cl
       Ptoto=Ptotl
       etoto=etotl
       vdotBo=vdotBl
       emago=emagl
    else if(SAL>0d0)then
       icase=2
       ro=rstarl
       uo=ustar
       vo=vstarl
       wo=wstarl
       Bo=Bstarl
       Co=Cstarl
       Ptoto=Ptotstar
       etoto=etotstarl
       vdotBo=vdotBstarl
       emago=half*(A*A+Bo*Bo+Co*Co)
    else if(ustar>0d0)then
       icase=3
       ro=rstarl
       uo=ustar
       vo=vstarstar
       wo=wstarstar
       Bo=Bstarstar
       Co=Cstarstar
       Ptoto=Ptotstar
       etoto=etotstarstarl
       vdotBo=vdotBstarstar
       emago=half*(A*A+Bo*Bo+Co*Co)
    else if(SAR>0d0)then
       icase=4
       ro=rstarr
       uo=ustar
       vo=vstarstar
       wo=wstarstar
       Bo=Bstarstar
       Co=Cstarstar
       Ptoto=Ptotstar
       etoto=etotstarstarr
       vdotBo=vdotBstarstar
       emago=half*(A*A+Bo*Bo+Co*Co)
    else if (SR>0d0)then
       icase=5
       ro=rstarr
       uo=ustar
       vo=vstarr
       wo=wstarr
       Bo=Bstarr
       Co=Cstarr
       Ptoto=Ptotstar
       etoto=etotstarr
       vdotBo=vdotBstarr
       emago=half*(A*A+Bo*Bo+Co*Co)
    else
       icase=6
       ro=rr
       uo=ur
       vo=vr
       wo=wr
       Bo=Br
       Co=Cr
       Ptoto=Ptotr
       etoto=etotr
       vdotBo=vdotBr
       emago=emagr
    end if
    !---------------------------------------------------------------------------
    ! Add a diffusive flux in the down grad(ln(rho)) direction, which is only
    ! nonzero when the smallest of the left/right density is within a factor 10
    ! from smallr, and which has at most magnitude of 0.25 times the largest of
    ! the left/right fast mode speeds.
    !---------------------------------------------------------------------------
    ustar = alog(rl/rr)*max(0.0, alog10((10.*smallr)/min(rl,rr)))
    ustar = 0.25*max(cfastr,cfastl)*min(1.0, max(-1.0, ustar))
    uo = uo + ustar
    !
    ! Compute the Godunov flux
    flux(i,r1,j1) = ro*uo
    flux(i,e1,j1) = ro*uo*uo+Ptoto+emago*(lim-1d0)-A*A*lim
    flux(i,v1,j1) = (etoto+Ptoto)*uo-A*vdotBo
    flux(i,v2,j1) = ro*uo*vo-A*Bo*lim
    flux(i,v3,j1) = ro*uo*wo-A*Co*lim
  end do
end subroutine hlld_ramses_mhd

!===============================================================================
!> Produce the array that may be used to loop over calls to the Riemann-2D
!> solver.  The input arrays should hold:
!>
!>    v: The primitive variables
!>  dv1: The 1st perpendicular increment, ds/2*dvdy in the x-case
!>  dv2: The 2nd perpendicular increment, ds/2*dvdz in the x-case
!>
!> i index     par    perp1   perp2
!> --------------------------------
!>   1:         x       y       z
!>   2:         y       z       x
!>   2:         z       x       y
!>
!> j index:  perp1  perp2   iy     iz
!> ----------------------------------
!>   1: SW   -dy/2  -dz/2   iy     iz
!>   2: NW   -dy/2  +dz/2   iy     iz-1
!>   3: NE   +dy/2  +dz/2   iy-1   iz-1
!>   4: SE   +dy/2  -dz/2   iy-1   iz
!===============================================================================
SUBROUTINE solve_2d (l, u, q, dq, e, a)
#ifdef OFFLOAD
  !$omp declare target
#endif
  integer:: l(4), u(4)
  real, contiguous, dimension(:,:,:,:)     :: q     ! cell center
  real, contiguous, dimension(:,:,:,:)     :: e     ! cell center
  real, contiguous, dimension(:,:,:,:,:)   :: dq    ! ds/2*dvdy
  real, contiguous, dimension(:,:,:)       :: a     ! edge values
  intent(in):: q, dq
  !.............................................................................
  integer:: i1, i2, i3, i4, j1, j2, j3
  integer:: i21, i31, i41
  integer:: i22, i32, i42
  integer:: i23, i33, i43
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%solve_2d', itimer=itimer)
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    do j3=1,3                             ! axis 3
      j1 = 1 + mod(j3,3)                  ! axis 1
      j2 = 1 + mod(j1,3)                  ! axis 2
      !-------------------------------------------------------------------------
      ! For cell-centered variables we get to the edges by down-staggering in
      ! the two directions perpendicular to the E-direction being computed.
      !-------------------------------------------------------------------------
      do i4=l(4),u(4)
        select case (j3)
        case (1)
          i21 = 0; i22 = 0; i23 = 1         ! LR-edge
          i31 = 0; i32 = 1; i33 = 1         ! RR-edge
          i41 = 0; i42 = 1; i43 = 0         ! RL-edge
        case (2)
          i21 = 1; i22 = 0; i23 = 0         ! LR-edge
          i31 = 1; i32 = 0; i33 = 1         ! RR-edge
          i41 = 0; i42 = 0; i43 = 1         ! RL-edge
        case (3)
          i21 = 0; i22 = 1; i23 = 0         ! LR-edge
          i31 = 1; i32 = 1; i33 = 0         ! RR-edge
          i41 = 1; i42 = 0; i43 = 0         ! RL-edge
        end select
        !-----------------------------------------------------------------------
        ! At this point all primitive values in the q() array need to be the
        ! ones predicted at t+dt/2, not the original ones. The values that come
        ! out of this loop are good to uo = ub-2.
        !-----------------------------------------------------------------------
        !$omp simd
        do i1=l(1),u(1)
          !---------------------------------------------------------------------
          ! When doing Ez, this case has no index shifts
          !---------------------------------------------------------------------
          a(i1,i4,1) =  &                                ! RR-edge
                        q(i1    ,i4   ,i2    ,i3    )  & ! cell center
                     - dq(i1    ,i4,j1,i2    ,i3    )  & ! dx/2 increment
                     - dq(i1    ,i4,j2,i2    ,i3    )    ! dy/2 increment
          !---------------------------------------------------------------------
          ! When doing Ez, this case an index shift down in y; i22=1
          !---------------------------------------------------------------------
          a(i1,i4,2) =  &                                ! RL-edge
                        q(i1-i21,i4   ,i2-i22,i3-i23)  & ! cell center
                     - dq(i1-i21,i4,j1,i2-i22,i3-i23)  & ! dx/2 increment
                     + dq(i1-i21,i4,j2,i2-i22,i3-i23)    ! dy/2 increment
          !---------------------------------------------------------------------
          ! When doing Ez, this case index shifts down in x & y; i31=i32=1
          !---------------------------------------------------------------------
          a(i1,i4,3) =  &                                ! LL-edge
                        q(i1-i31,i4   ,i2-i32,i3-i33)  & ! cell center
                     + dq(i1-i31,i4,j1,i2-i32,i3-i33)  & ! dy/2 increment
                     + dq(i1-i31,i4,j2,i2-i32,i3-i33)    ! dz/2 increment
          !---------------------------------------------------------------------
          ! When doing Ez, this case an index shift down in x; i41=1
          !---------------------------------------------------------------------
          a(i1,i4,4) =  &                                ! LR-edge
                        q(i1-i41,i4   ,i2-i42,i3-i43)  & ! cell center
                     + dq(i1-i41,i4,j1,i2-i42,i3-i43)  & ! dx/2 increment
                     - dq(i1-i41,i4,j2,i2-i42,i3-i43)    ! dy/2 increment
        end do
      end do
      !-------------------------------------------------------------------------
      ! Get denisty and pressure values from the log values
      !-------------------------------------------------------------------------
      !$omp simd
      do i1=l(1),u(1)
        a(i1,r1,1) = exp(a(i1,r1,1))
        a(i1,r1,2) = exp(a(i1,r1,2))
        a(i1,r1,3) = exp(a(i1,r1,3))
        a(i1,r1,4) = exp(a(i1,r1,4))
        a(i1,e1,1) = exp(a(i1,e1,1))
        a(i1,e1,2) = exp(a(i1,e1,2))
        a(i1,e1,3) = exp(a(i1,e1,3))
        a(i1,e1,4) = exp(a(i1,e1,4))
      end do
      call hlld_2d (l,u,j3,a(:,:,1),a(:,:,2),a(:,:,3),a(:,:,4), e(:,i2,i3,j3))
    end do
  end do
  end do
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE solve_2d

!===============================================================================
!> 2-D HLLD Riemann solver, with four edge values already properly selected,
!> but with variable indices in the original order; reordered below
!===============================================================================
SUBROUTINE hlld_2d (lb, ub, j, qRR, qRL, qLL, qLR, emf)
#ifdef OFFLOAD
  !$omp declare target
#endif
  integer, intent(in):: lb(4), ub(4), j
  real, contiguous, dimension(:,:):: qLL,qRL,qLR,qRR
  real, contiguous, dimension(:):: emf
  intent(in):: qLL, qRL, qLR, qRR
  !.............................................................................
  ! local variables
  real:: e,ell,erl,elr,err,s,sl,sr,sb,st,sal,sar,sat,sab
  real:: cfastLLx,cfastRLx,cfastLRx,cfastRRx,cfastLLy,cfastRLy,cfastLRy,cfastRRy
  real:: calfvenR,calfvenL,calfvenT,calfvenB
  real:: vLLx,vRLx,vLRx,vRRx,vLLy,vRLy,vLRy,vRRy
  real:: PtotLL,PtotLR,PtotRL,PtotRR,rcLLx,rcLRx,rcRLx,rcRRx,rcLLy,rcLRy,rcRLy,rcRRy
  real:: ustar,vstar,rstarLLx,rstarLRx,rstarRLx,rstarRRx,rstarLLy,rstarLRy,rstarRLy,rstarRRy
  real:: rstarLL,rstarLR,rstarRL,rstarRR,AstarLL,AstarLR,AstarRL,AstarRR,BstarLL,BstarLR,BstarRL,BstarRR
  real:: EstarLLx,EstarLRx,EstarRLx,EstarRRx,EstarLLy,EstarLRy,EstarRLy,EstarRRy,EstarLL,EstarLR,EstarRL,EstarRR
  real:: AstarT,AstarB,BstarR,BstarL
  !
  real:: rLL,rLR,rRL,rRR,pLL,pLR,pRL,pRR,uLL,uLR,uRL,uRR,vLL,vLR,vRL,vRR, &
         ALL,ALR,ARL,ARR,BLL,BLR,BRL,BRR,CLL,CLR,CRL,CRR
  real:: qtmp1,qtmp2,qtmp4,q,qtmp6,qtmp8
  real:: b2, c2, d2, d1, gam
  integer:: i, j1, j2, j3, v1, v2, v3, n1, n2, n3
  !-----------------------------------------------------------------------------
  ! j1 is the E-component direction, with edge values qLL,qLR,qRR,qRL
  ! j2,j3 are the two perpendicular directions.  v1-v3 and n1-n3 are the
  ! corresponding velocity and magnetic field components.
  !-----------------------------------------------------------------------------
  j3 = j
  j1 = 1 + mod(j3,3)
  j2 = 1 + mod(j1,3)
  v1 = u1+j1-1
  v2 = u1+j2-1
  v3 = u1+j3-1
  n1 = m1+j1-1
  n2 = m1+j2-1
  n3 = m1+j3-1
  gam = gamma
  !$omp simd
  do i=lb(1),ub(1)
    ! from the code below it is clear that in the simple limit the E-component
    ! is evaluated as E = u*B - v*A 
    rLL=qLL(i,r1); pLL=qLL(i,e1)
    uLL=qLL(i,v1); vLL=qLL(i,v2)
    ALL=qLL(i,n1); BLL=qLL(i,n2); CLL=qLL(i,n3)
    !
    rLR=qLR(i,r1); pLR=qLR(i,e1)
    uLR=qLR(i,v1); vLR=qLR(i,v2) 
    ALR=qLR(i,n1); BLR=qLR(i,n2); CLR=qLR(i,n3)
    !
    rRL=qRL(i,r1); pRL=qRL(i,e1)
    uRL=qRL(i,v1); vRL=qRL(i,v2)
    ARL=qRL(i,n1); BRL=qRL(i,n2); CRL=qRL(i,n3)
    !
    rRR=qRR(i,r1); pRR=qRR(i,e1)
    uRR=qRR(i,v1); vRR=qRR(i,v2)
    ARR=qRR(i,n1); BRR=qRR(i,n2); CRR=qRR(i,n3)
    !---------------------------------------------------------------------------
    ! Compute 4 fast magnetosonic velocities relative to x direction
    !---------------------------------------------------------------------------
    b2 = ALL**2+BLL**2+CLL**2
    d1 = 1./rLL
    c2 = gam*pLL*d1
    d2 = half*(b2*d1+c2)
    cfastLLx=sqrt(d2 + sqrt(max(d2*d2-c2*ALL**2*d1,zero)))
    !
    b2 = ARL**2+BRL**2+CRL**2
    d1 = 1./rRL
    c2 = gam*pRL*d1
    d2 = half*(b2*d1+c2)
    cfastRLx=sqrt(d2 + sqrt(max(d2*d2-c2*ARL**2*d1,zero)))
    !
    b2 = ALR**2+BLR**2+CLR**2
    d1 = 1./rLR
    c2 = gam*pLR*d1
    d2 = half*(b2*d1+c2)
    cfastLRx=sqrt(d2 + sqrt(max(d2*d2-c2*ALR**2*d1,zero)))
    !
    b2 = ARR**2+BRR**2+CRR**2
    d1 = 1./rRR
    c2 = gam*pRR*d1
    d2 = half*(b2*d1+c2)
    cfastRRx=sqrt(d2 + sqrt(max(d2*d2-c2*ARR**2*d1,zero)))
    !
    !---------------------------------------------------------------------------
    ! Compute 4 fast magnetosonic velocity relative to y direction
    !---------------------------------------------------------------------------
    b2 = ALL**2+BLL**2+CLL**2
    d1 = 1./rLL
    c2 = gam*pLL*d1
    d2 = half*(b2*d1+c2)
    cfastLLy=sqrt(d2 + sqrt(max(d2*d2-c2*BLL**2*d1,zero)))
    !
    b2 = ALR**2+BLR**2+CLR**2
    d1 = 1./rLR
    c2 = gam*pLR*d1
    d2 = half*(b2*d1+c2)
    cfastLRy=sqrt(d2 + sqrt(max(d2*d2-c2*BLR**2*d1,zero)))
    !
    b2 = ARL**2+BRL**2+CRL**2
    d1 = 1./rRL
    c2 = gam*pRL*d1
    d2 = half*(b2*d1+c2)
    cfastRLy=sqrt(d2 + sqrt(max(d2*d2-c2*BRL**2*d1,zero)))
    !
    b2 = ARR**2+BRR**2+CRR**2
    d1 = 1./rRR
    c2 = gam*pRR*d1
    d2 = half*(b2*d1+c2)
    cfastRRy=sqrt(d2 + sqrt(max(d2*d2-c2*BRR**2*d1,zero)))
    !
    !
    SL=min(uLL,uLR,uRL,uRR)-max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
    SR=max(uLL,uLR,uRL,uRR)+max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
    SB=min(vLL,vLR,vRL,vRR)-max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
    ST=max(vLL,vLR,vRL,vRR)+max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
    !
    if(SB > 0.0) then
       if(SL > 0.0) then
          ELL=uLL*BLL-vLL*ALL
          E=ELL
       else if(SR < 0.0) then
          ERL=uRL*BRL-vRL*ARL
          E=ERL
       else
          PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
          PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
          PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
          PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
          rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
          rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
          rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
          rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
    !
  ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR)) &
        /(rcLLx+rcLRx+rcRLx+rcRRx)
  vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR)) &
        /(rcLLy+rcLRy+rcRLy+rcRRy)
    !
          rstarLLx=rLL*(SL-uLL)/(SL-ustar)
          rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
          AstarLL=ALL*(SB-vLL)/(SB-vstar)
    !
          rstarLRx=rLR*(SL-uLR)/(SL-ustar)
          rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
          AstarLR=ALR*(ST-vLR)/(ST-vstar)
    !
          rstarRLx=rRL*(SR-uRL)/(SR-ustar)
          rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
          AstarRL=ARL*(SB-vRL)/(SB-vstar)
    !
          rstarRRx=rRR*(SR-uRR)/(SR-ustar)
          rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
          AstarRR=ARR*(ST-vRR)/(ST-vstar)
    !
          calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                       abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
          calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                       abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
    !
          SAL=min(ustar-calfvenL,zero)
          SAR=max(ustar+calfvenR,zero)
    !
          BstarLL=BLL*(SL-uLL)/(SL-ustar)
          BstarRL=BRL*(SR-uRL)/(SR-ustar)
    !
          EstarLLx=ustar*BstarLL-vLL*ALL
          EstarRLx=ustar*BstarRL-vRL*ARL
    !
          E=(SAR*EstarLLx-SAL*EstarRLx+SAR*SAL*(BRL-BLL))/(SAR-SAL)
       endif
    else if (ST < 0.0) then
      if(SL > 0.0)then
          ELR=uLR*BLR-vLR*ALR
    !
          E=ELR
      else if(SR < 0.0) then
          ERR=uRR*BRR-vRR*ARR
    !
          E=ERR
      else
          PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
          PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
          PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
          PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
    !
          rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
          rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
          rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
          rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
    !
  ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR)) &
        /(rcLLx+rcLRx+rcRLx+rcRRx)
  vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR)) &
        /(rcLLy+rcLRy+rcRLy+rcRRy)
    !
          rstarLLx=rLL*(SL-uLL)/(SL-ustar)
          rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
          AstarLL=ALL*(SB-vLL)/(SB-vstar)
    !
          rstarLRx=rLR*(SL-uLR)/(SL-ustar)
          rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
          AstarLR=ALR*(ST-vLR)/(ST-vstar)
    !
          rstarRLx=rRL*(SR-uRL)/(SR-ustar)
          rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
          AstarRL=ARL*(SB-vRL)/(SB-vstar)
    !
          rstarRRx=rRR*(SR-uRR)/(SR-ustar)
          rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
          AstarRR=ARR*(ST-vRR)/(ST-vstar)
    !
          calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                       abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
          calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                       abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
    !
          SAL=min(ustar-calfvenL,zero)
          SAR=max(ustar+calfvenR,zero)
    !
          BstarLR=BLR*(SL-uLR)/(SL-ustar)
          BstarRR=BRR*(SR-uRR)/(SR-ustar)
    !
          EstarLRx=ustar*BstarLR-vLR*ALR
          EstarRRx=ustar*BstarRR-vRR*ARR
    !
          E=(SAR*EstarLRx-SAL*EstarRRx+SAR*SAL*(BRR-BLR))/(SAR-SAL)
       endif
    else if(SL > 0.0)then
       PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
       PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
       PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
       PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
    !
       rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
       rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
       rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
       rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
    !
  ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR)) &
      /(rcLLx+rcLRx+rcRLx+rcRRx)
  vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR)) &
      /(rcLLy+rcLRy+rcRLy+rcRRy)
    !
       rstarLLy=rLL*(SB-vLL)/(SB-vstar)
       rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
       BstarLL=BLL*(SL-uLL)/(SL-ustar)
    !
       rstarLRy=rLR*(ST-vLR)/(ST-vstar)
       rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
       BstarLR=BLR*(SL-uLR)/(SL-ustar)
    !
       rstarRLy=rRL*(SB-vRL)/(SB-vstar)
       rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
       BstarRL=BRL*(SR-uRL)/(SR-ustar)
    !
       rstarRRy=rRR*(ST-vRR)/(ST-vstar)
       rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
       BstarRR=BRR*(SR-uRR)/(SR-ustar)
    !
       calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                    abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
       calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                    abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
    !
       SAB=min(vstar-calfvenB,zero)
       SAT=max(vstar+calfvenT,zero)
    !
       AstarLL=ALL*(SB-vLL)/(SB-vstar)
       AstarLR=ALR*(ST-vLR)/(ST-vstar)
    !
       EstarLLy=uLL*BLL-vstar*AstarLL
       EstarLRy=uLR*BLR-vstar*AstarLR
    !
       E=(SAT*EstarLLy-SAB*EstarLRy-SAT*SAB*(ALR-ALL))/(SAT-SAB)
    else if(SR < 0.0)then
       PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
       PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
       PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
       PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
    !
       rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
       rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
       rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
       rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
    !
  ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR)) &
      /(rcLLx+rcLRx+rcRLx+rcRRx)
  vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR)) &
      /(rcLLy+rcLRy+rcRLy+rcRRy)
    !
       rstarLLy=rLL*(SB-vLL)/(SB-vstar)
       rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
       BstarLL=BLL*(SL-uLL)/(SL-ustar)
    !
       rstarLRy=rLR*(ST-vLR)/(ST-vstar)
       rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
       BstarLR=BLR*(SL-uLR)/(SL-ustar)
    !
       rstarRLy=rRL*(SB-vRL)/(SB-vstar)
       rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
       BstarRL=BRL*(SR-uRL)/(SR-ustar)
    !
       rstarRRy=rRR*(ST-vRR)/(ST-vstar)
       rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
       BstarRR=BRR*(SR-uRR)/(SR-ustar)
    !
       calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                    abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
       calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                    abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
    !
       SAB=min(vstar-calfvenB,zero)
       SAT=max(vstar+calfvenT,zero)
    !
       AstarRL=ARL*(SB-vRL)/(SB-vstar)
       AstarRR=ARR*(ST-vRR)/(ST-vstar)
    !
       EstarRLy=uRL*BRL-vstar*AstarRL
       EstarRRy=uRR*BRR-vstar*AstarRR
    !
       E=(SAT*EstarRLy-SAB*EstarRRy-SAT*SAB*(ARR-ARL))/(SAT-SAB)
    else
       PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
       PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
       PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
       PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
    !
       rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
       rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
       rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
       rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
    !
       ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))&
            /(rcLLx+rcLRx+rcRLx+rcRRx)
       vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))&
       /(rcLLy+rcLRy+rcRLy+rcRRy)
    !
       rstarLLx=rLL*(SL-uLL)/(SL-ustar); BstarLL=BLL*(SL-uLL)/(SL-ustar)
       rstarLLy=rLL*(SB-vLL)/(SB-vstar); AstarLL=ALL*(SB-vLL)/(SB-vstar)
       rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
    !
       rstarLRx=rLR*(SL-uLR)/(SL-ustar); BstarLR=BLR*(SL-uLR)/(SL-ustar)
       rstarLRy=rLR*(ST-vLR)/(ST-vstar); AstarLR=ALR*(ST-vLR)/(ST-vstar)
       rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
    !
       rstarRLx=rRL*(SR-uRL)/(SR-ustar); BstarRL=BRL*(SR-uRL)/(SR-ustar)
       rstarRLy=rRL*(SB-vRL)/(SB-vstar); AstarRL=ARL*(SB-vRL)/(SB-vstar)
       rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
    !
       rstarRRx=rRR*(SR-uRR)/(SR-ustar); BstarRR=BRR*(SR-uRR)/(SR-ustar)
       rstarRRy=rRR*(ST-vRR)/(ST-vstar); AstarRR=ARR*(ST-vRR)/(ST-vstar)
       rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
    !
       calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                    abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
       calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                    abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
       calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                    abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
       calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                    abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
       SAL=min(ustar-calfvenL,zero)
       SAR=max(ustar+calfvenR,zero)
       SAB=min(vstar-calfvenB,zero)
       SAT=max(vstar+calfvenT,zero)
       EstarLL =ustar*BstarLL-vstar*AstarLL
       EstarRL =ustar*BstarRL-vstar*AstarRL
       EstarLR =ustar*BstarLR-vstar*AstarLR
       EstarRR =ustar*BstarRR-vstar*AstarRR
    !
       AstarT=(SAR*AstarRR-SAL*AstarLR)/(SAR-SAL)
       AstarB=(SAR*AstarRL-SAL*AstarLL)/(SAR-SAL)
       BstarR=(SAT*BstarRR-SAB*BstarRL)/(SAT-SAB)
       BstarL=(SAT*BstarLR-SAB*BstarLL)/(SAT-SAB)
    !
       E = (SAL*SAB*EstarRR-SAL*SAT*EstarRL &
           -SAR*SAB*EstarLR+SAR*SAT*EstarLL)/(SAR-SAL)/(SAT-SAB) &
           -SAT*SAB/(SAT-SAB)*(AstarT-AstarB) &
           +SAR*SAL/(SAR-SAL)*(BstarR-BstarL)
    endif
    emf(i) = -E
  end do
END SUBROUTINE hlld_2d

!===============================================================================
!> Use the divergence of the HD fluxes and the curl of the EMFs to update the
!> conservative variables.
!===============================================================================
  SUBROUTINE div_flux (self, l, u, uold, unew, flux, emf, dtds)
#ifdef OFFLOAD
  !$omp declare target
#endif
  class(mhd_t):: self
  integer:: l(4), u(4)
  real, contiguous:: uold(:,:,:,:), unew(:,:,:,:), flux(:,:,:,:,:), emf(:,:,:,:)
  real:: dtds
  !.............................................................................
  real:: d, e, px, py, pz, bx, by, bz
  real:: f, dmin0, dmin1, g0, g1, ek0, ek1, em0, em1, eth, eth_min, d_min, p_min
  real:: ux, uy, uz, d0, d1, eth0, eth1, dtds1
  integer:: i1, i2, i3, i4, n_smalld, n_smallp, n_iter, unit
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
  n_iter = 0
  do while (f < 1.0-courant_d .and. n_iter<10)
    self%dtime = f*self%dtime
    dtds1      = f*dtds1
    if (verbose > 0 .and. n_iter==0) then
      write (stdout,'(a,3i6,f7.3,1p,3e12.3)') &
        'div_flux: id,ist,f,d_min,p_min,dtds_f=', &
        self%id, n_iter, self%istep, f, d_min, p_min, dtds1/dtds
    end if
    call divfl
    if (verbose > 0) then
      write (stdout,'(a,3i6,f7.3,1p,3e12.3)') &
        'div_flux: id,ist,f,d_min,p_min,dtds_f=', &
        self%id, n_iter, self%istep, f, d_min, p_min, dtds1/dtds
    end if
    n_iter = n_iter+1
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
#ifdef OFFLOAD
  !$omp parallel do collapse(2) schedule(static,1)
#endif
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !$omp simd
    do i1=l(1),u(1)
      unew(i1,i2,i3,m1) = &
      uold(i1,i2,i3,m1) - dtds1*(emf(i1,i2+1,i3,3) - emf(i1,i2,i3,3) &
                               - emf(i1,i2,i3+1,2) + emf(i1,i2,i3,2))
      unew(i1,i2,i3,m2) = &
      uold(i1,i2,i3,m2) - dtds1*(emf(i1,i2,i3+1,1) - emf(i1,i2,i3,1) &
                               - emf(i1+1,i2,i3,3) + emf(i1,i2,i3,3))
      unew(i1,i2,i3,m3) = &
      uold(i1,i2,i3,m3) - dtds1*(emf(i1+1,i2,i3,2) - emf(i1,i2,i3,2) &
                               - emf(i1,i2+1,i3,1) + emf(i1,i2,i3,1))
    end do
  end do
  end do
  n_smalld = 0
  n_smallp = 0
  f = 1.0
  d_min = 1e9
  p_min = 1e9
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i4=1,5
        !$omp simd
        do i1=l(1),u(1)
          unew(i1,i2,i3,i4) = &
          uold(i1,i2,i3,i4) + dtds1*( &
            flux(i1,i4,1,i2,i3) - flux(i1+1,i4,1,i2  ,i3  ) &
          + flux(i1,i4,2,i2,i3) - flux(i1  ,i4,2,i2+1,i3  ) &
          + flux(i1,i4,3,i2,i3) - flux(i1  ,i4,3,i2  ,i3+1))
        end do
      end do
      do i1=l(1),u(1)
        d0 = uold(i1,i2,i3,1)
        ux = uold(i1,i2,i3,2)/d0
        uy = uold(i1,i2,i3,3)/d0
        uz = uold(i1,i2,i3,4)/d0
        ek0 = 0.5*d0*(ux**2+uy**2+uz**2)
        bx = 0.5*(uold(i1,i2,i3,6) + uold(i1+1,i2,i3,6))
        by = 0.5*(uold(i1,i2,i3,7) + uold(i1,i2+1,i3,7))
        bz = 0.5*(uold(i1,i2,i3,8) + uold(i1,i2,i3+1,8))
        em0 = 0.5*(bx**2+by**2+bz**2)
        eth0 = uold(i1,i2,i3,5)
        d1 = unew(i1,i2,i3,1)
        d_min = min(d_min,d1)
        n_smalld = merge(n_smalld+1,n_smalld,d1<smallr)
        d1 = max(d1,smallr)
        unew(i1,i2,i3,1) = d1
        ux = unew(i1,i2,i3,2)/d1
        uy = unew(i1,i2,i3,3)/d1
        uz = unew(i1,i2,i3,4)/d1
        ek1 = 0.5*d1*(ux**2+uy**2+uz**2)
        bx = 0.5*(unew(i1,i2,i3,6) + unew(i1+1,i2,i3,6))
        by = 0.5*(unew(i1,i2,i3,7) + unew(i1,i2+1,i3,7))
        bz = 0.5*(unew(i1,i2,i3,8) + unew(i1,i2,i3+1,8))
        em1 = 0.5*(bx**2+by**2+bz**2)
        eth1 = unew(i1,i2,i3,5) + (ek0-ek1) + (em0-em1)
        eth_min = smallc**2/(g1*g0)*d1**g0
        p_min = min(p_min,eth_min*g1*g0)
        n_smallp = merge(n_smallp+1,n_smallp,eth1<eth_min)
        eth1 = max(eth1,eth_min)
        unew(i1,i2,i3,5) = eth1
        f = min(f,d1/d0,eth1/eth0)
      end do
    end do
  end do
  f = exp(-(1.0-f))
end subroutine divfl
END SUBROUTINE div_flux

!===============================================================================
!> Add the effect of source terms.  Since we keep only the thermal energy in
!> the time slices, the only effect of force terms is to change the linear
!> momementa, and the only effect of the heating/cooling is to change the
!> thermal energy.
!===============================================================================
SUBROUTINE sources (l, u, unew, dt)
  integer :: l(3), u(3)
  real    :: unew(:,:,:,:)
  real    :: dt
  !.............................................................................
  integer :: i1, i2, i3
  real    :: d
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('mhd_t%sources', itimer=itimer)
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
